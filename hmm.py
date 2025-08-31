#!/usr/bin/env python3
import argparse, os, re, subprocess, sys, shlex, tempfile, shutil
from collections import defaultdict, OrderedDict

import numpy as np
from pyfaidx import Fasta

# ----------------- helpers -----------------

def run(cmd, check=True, capture=False):
    print("[cmd]", cmd)
    if capture:
        p = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        if check and p.returncode != 0:
            sys.stderr.write(p.stdout); sys.stderr.write(p.stderr)
            raise RuntimeError(f"Command failed: {cmd}")
        return p.stdout
    else:
        p = subprocess.run(cmd, shell=True)
        if check and p.returncode != 0:
            raise RuntimeError(f"Command failed: {cmd}")
        return ""

def parse_region(s):
    # CHR:START-END (1-based, inclusive)
    m = re.match(r"^\s*([^:]+):\s*([0-9,._]+)\s*[-:]\s*([0-9,._]+)\s*$", s)
    if not m:
        raise ValueError(f"Bad --region: {s}")
    chrom = m.group(1)
    start = int(re.sub(r"[,_\.]", "", m.group(2)))
    end   = int(re.sub(r"[,_\.]", "", m.group(3)))
    if start < 1 or end < start:
        raise ValueError(f"Bad region bounds: {start}-{end}")
    return chrom, start, end

def tersect_view(tsi, expr, region, extra=""):
    # returns list of VCF-like rows as strings
    cmd = f"tersect view {shlex.quote(tsi)} \"{expr}\" {shlex.quote(region)} {extra}"
    out = run(cmd, capture=True)
    lines = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
    return lines

def convert_tree_to_species(tree_newick, sample_to_species):
    """Convert a tree with sample names to species names"""
    converted_tree = tree_newick
    
    # Replace each sample name with its corresponding species name
    for sample, species in sample_to_species.items():
        # Use word boundaries to avoid partial matches
        import re
        pattern = r'\b' + re.escape(sample) + r'\b'
        converted_tree = re.sub(pattern, species, converted_tree)
    
    return converted_tree

def parse_vcf_like(lines):
    # Expect at least CHROM, POS, REF, ALT columns (TSI prints VCF-like)
    recs = []
    for ln in lines:
        cols = re.split(r"\s+", ln.strip())
        if len(cols) < 5:  # CHROM POS ID REF ALT ...
            continue
        chrom, pos, _id, ref, alt = cols[:5]
        try:
            pos_i = int(pos)
        except ValueError:
            continue
        recs.append((chrom, pos_i, ref, alt))
    return recs

def load_panel(tsi, samples, region, biallelic_only=True):
    # Process each sample individually since union doesn't work
    all_recs = []
    
    for sample in samples:
        print(f"[info] Processing sample: {sample}")
        lines = tersect_view(tsi, f"'{sample}'", region)
        sample_recs = parse_vcf_like(lines)
        all_recs.extend(sample_recs)
    
    # Now merge all records by position
    by_pos = defaultdict(list)
    for chrom, pos, ref, alt in all_recs:
        if len(ref) == 1 and len(alt) == 1 and ref in "ACGTacgt" and alt in "ACGTacgt":
            by_pos[(chrom, pos)].append((ref.upper(), alt.upper()))
    
    panel = []
    for (chrom, pos), alleles in sorted(by_pos.items()):
        alts = {a for (_, a) in alleles}
        refs = {r for (r, _) in alleles}
        if biallelic_only and (len(alts) != 1 or len(refs) != 1):
            continue  # skip multiallelic or inconsistent REF
        # pick the majority ALT / single ALT
        alt = next(iter(alts))
        ref = next(iter(refs)) if refs else None
        panel.append((chrom, pos, ref, alt))
    
    return panel  # list of (chrom, pos, ref, alt)

def fasta_ref_base(fa, chrom, pos):
    try:
        base = fa[chrom][pos-1].seq.upper()
    except KeyError:
        raise SystemExit(f"Chromosome not in FASTA: {chrom}")
    except Exception as e:
        raise SystemExit(f"FASTA fetch failed at {chrom}:{pos}: {e}")
    if base not in "ACGT":
        return "N"
    return base

def build_sample_alt_positions(tsi, sample, region):
    # set of positions where sample has ALT allele(s)
    print(f"[info] Building ALT positions for sample: {sample}")
    lines = tersect_view(tsi, f"'{sample}'", region)
    recs = parse_vcf_like(lines)
    return {(chrom, pos): alt.upper() for chrom, pos, ref, alt in recs}

def infer_species_name(sample, tokens=2):
    # Heuristic: join first <tokens> underscore-separated pieces
    parts = sample.split("_")
    if len(parts) >= tokens:
        return "_".join(parts[:tokens])
    return parts[0]

def write_fasta(path, sequences):
    # sequences: OrderedDict {sample -> seq}
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

def write_nexus_hmm(out_path, sequences, network_newick, allelemap, model="gtr",
                    iterations=20, threads=4, runs=1, outdir="phylonet_hmm_out"):
    samples = list(sequences.keys())
    Lset = {len(sequences[s]) for s in samples}
    if len(Lset) != 1:
        raise SystemExit(f"Sequences must be same length; saw lengths {sorted(Lset)}")
    nchar = Lset.pop()
    ntax = len(samples)

    # allele map format: <Sp1:s1,s2; Sp2:s3>
    # Fixed: use commas between samples of same species, filter empty species
    parts = []
    for sp, sample_list in allelemap.items():
        if sample_list:  # Only include species that have samples
            parts.append(f"{sp}:{','.join(sample_list)}")  # Use commas, not spaces
    amap = "<" + "; ".join(parts) + ">"

    model_flag = "-gtr" if model.lower() == "gtr" else "-jc"

    with open(out_path, "w") as nex:
        nex.write("#NEXUS\n\n")
        nex.write("BEGIN NETWORKS;\n")
        nex.write(f"  Network net = {network_newick.rstrip(';')};\n")
        nex.write("END;\n\n")
        nex.write("BEGIN DATA;\n")
        nex.write(f"  dimensions ntax={ntax} nchar={nchar};\n")
        nex.write('  format datatype=dna symbols="ACGT" missing=? gap=-;\n')
        nex.write("  matrix\n")
        for s in samples:
            nex.write(f"    {s:<20} {sequences[s]}\n")
        nex.write("  ;\nEND;\n\n")
        nex.write("BEGIN PHYLONET;\n")
        nex.write(
            f"  HmmCommand net -allelemap {amap} "
            f"{model_flag} -iterations {20} -threads {14} "
            f"-numberofruns {1} -outputdirectory \"{outdir}\";\n"
        )
        nex.write("END;\n")
def ensure_iqtree(aln_fa, outdir, threads=4):
    # use iqtree2 if available, else iqtree
    exe = None
    for cand in ("iqtree2", "iqtree"):
        if shutil.which(cand):
            exe = cand
            break
    if exe is None:
        raise SystemExit("IQ-TREE not found. Install conda package 'iqtree' (or 'iqtree2').")
    
    # Remove -B 0 since it's causing issues, just do ML tree without bootstrap
    cmd = f"{exe} -s {shlex.quote(aln_fa)} -m GTR+G -T {threads} -redo -quiet"
    run(cmd)
    treefile = aln_fa + ".treefile"
    if not os.path.exists(treefile):
        raise SystemExit("IQ-TREE did not produce .treefile")
    
    with open(treefile, 'r') as f:
        return f.read().strip()

# ----------------- main -----------------

def main():
    ap = argparse.ArgumentParser(description="Build PhyloNet-HMM NEXUS from FASTA+TSI for a region.")
    ap.add_argument("--tsi", required=True, help="Path to Tersect .tsi")
    ap.add_argument("--fasta", required=True, help="Reference FASTA (indexed).")
    ap.add_argument("--region", required=True, help="Region CHR:START-END (1-based).")
    ap.add_argument("--samples", nargs="+", required=True, help="Sample names to include.")
    ap.add_argument("--species-map", help="TSV/CSV with columns: sample,species. If omitted, infer by name prefix.")
    ap.add_argument("--species-prefix-tokens", type=int, default=2, help="If inferring species, tokens to join (default 2).")
    ap.add_argument("--biallelic-only", action="store_true", default=True, help="Skip multiallelic sites (default).")
    ap.add_argument("--keep-multiallelic", dest="biallelic_only", action="store_false", help="Allow multiallelic (encodes as '?').")
    ap.add_argument("--model", choices=["gtr","jc"], default="gtr")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--iterations", type=int, default=300)
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--network", help="Extended Newick (.enewick). If omitted, infer ML tree with IQ-TREE.")
    ap.add_argument("--phylonet-cmd", default="phylonet", help='Command to run PhyloNet (e.g. "phylonet" or "java -jar PhyloNet.jar")')
    ap.add_argument("--outdir", default="phylonet_hmm_out")
    ap.add_argument("--prefix", default="hmm")
    ap.add_argument("--write-only", action="store_true", help="Write outputs but do not run PhyloNet.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Region + FASTA
    chrom, start, end = parse_region(args.region)
    fa = Fasta(args.fasta, as_raw=True, rebuild=False)

    # 2) Site panel (union of samples)
    panel = load_panel(args.tsi, args.samples, args.region, biallelic_only=args.biallelic_only)
    if not panel:
        sys.exit("No SNPs found in region for selected samples.")
    # fill missing REF from FASTA if needed; skip if mismatch and biallelic_only
    fixed_panel = []
    for (c, pos, ref, alt) in panel:
        if ref is None or ref not in "ACGT":
            ref = fasta_ref_base(fa, c, pos)
        if ref not in "ACGT" or alt not in "ACGT":
            continue
        fixed_panel.append((c, pos, ref, alt))
    panel = fixed_panel

    # 3) Per-sample ALT sets
    sample_alt = {}
    for s in args.samples:
        sample_alt[s] = build_sample_alt_positions(args.tsi, s, args.region)

    # 4) Build SNP alignment (one character per site in panel)
    # order panel by (chrom,pos)
    panel = sorted(panel, key=lambda x: (x[0], x[1]))
    positions = [(c,p) for (c,p,_,_) in panel]
    ref_arr = np.array([ref for (_,_,ref,_) in panel], dtype='<U1')  # Fixed numpy dtype
    alt_arr = np.array([alt for (_,_,_,alt) in panel], dtype='<U1')  # Fixed numpy dtype
    sequences = OrderedDict()
    for s in args.samples:
        arr = ref_arr.copy()
        alt_calls = sample_alt[s]
        # set ALT where sample has ALT at the same (chrom,pos); otherwise REF.
        # if multiallelic was allowed and sample ALT != panel ALT, write '?'
        for i,(c,p,ref,alt) in enumerate(panel):
            sa = alt_calls.get((c,p))
            if sa is None:
                # absent → either REF or missing; here we assume callable REF
                continue
            if sa == alt:
                arr[i] = alt
            else:
                arr[i] = "?"  # conflicting ALT at same pos
        sequences[s] = "".join(arr.tolist())

    # Write alignment FASTA (for transparency + tree inference)
    aln_fa = os.path.join(args.outdir, f"{args.prefix}.snp.fasta")
    write_fasta(aln_fa, sequences)
    print(f"[ok] Wrote SNP alignment FASTA: {aln_fa}")

    # 5) Species map → allele map
    species_to_samples = defaultdict(list)
    if args.species_map:
        # simple TSV/CSV sniff
        sep = "\t" if args.species_map.endswith((".tsv",".txt")) else ","
        with open(args.species_map) as fh:
            header = fh.readline().strip().split(sep)
            try:
                i_samp = header.index("sample"); i_spec = header.index("species")
            except ValueError:
                sys.exit("--species-map must have columns: sample,species")
            for ln in fh:
                cols = ln.rstrip("\n").split(sep)
                if len(cols) <= max(i_samp,i_spec): continue
                species_to_samples[cols[i_spec]].append(cols[i_samp])
    else:
        for s in args.samples:
            sp = infer_species_name(s, tokens=args.species_prefix_tokens)
            # Fix: if species name equals sample name, it means the sample doesn't follow
            # the expected naming pattern, so treat it as its own species
            if sp == s:
                print(f"[warning] Sample '{s}' doesn't follow expected naming pattern, treating as species '{s}'")
            species_to_samples[sp].append(s)

    # Filter out any species that don't have samples in our sequence data
    species_to_samples = {sp: samples for sp, samples in species_to_samples.items() 
                        if any(sample in sequences for sample in samples)}
    
    print(f"[debug] All samples in sequences: {list(sequences.keys())}")
    print(f"[debug] Species to samples mapping:")
    for species, sample_list in species_to_samples.items():
        print(f"  {species}: {sample_list}")

    # 6) Network/tree
    if args.network:
        with open(args.network, 'r') as f:
            network_newick = f.read().strip()
    else:
        # Build tree using one representative per species
        print("[info] Building species tree using representatives...")
        
        # Create species-representative sequences
        species_sequences = OrderedDict()
        for species, sample_list in species_to_samples.items():
            # Use first sample as representative for each species
            representative = sample_list[0]
            if representative in sequences:
                species_sequences[species] = sequences[representative]
                print(f"[info] Using {representative} as representative for {species}")
        
        print(f"[info] Found {len(species_sequences)} species: {list(species_sequences.keys())}")
        
        # Check if we have enough species for tree building
        if len(species_sequences) < 3:
            print(f"[warning] Only {len(species_sequences)} species found. PhyloNet-HMM needs at least 3 species.")
            if len(species_sequences) == 2:
                # Create a simple 2-species tree
                species_list = list(species_sequences.keys())
                network_newick = f"({species_list[0]},{species_list[1]});"
                print(f"[info] Using simple 2-species tree: {network_newick}")
            else:
                sys.exit("Need at least 2 species for PhyloNet-HMM analysis")
        else:
            # Write species alignment FASTA
            species_aln_fa = os.path.join(args.outdir, f"{args.prefix}.species.fasta")
            write_fasta(species_aln_fa, species_sequences)
            print(f"[ok] Wrote species alignment: {species_aln_fa}")
            
            # Build tree on species alignment
            network_newick = ensure_iqtree(species_aln_fa, args.outdir, threads=args.threads)

    print(f"[info] Using tree: {network_newick}")

    # 7) NEXUS for PhyloNet-HMM
    nexus_path = os.path.join(args.outdir, f"{args.prefix}.phylonet_hmm.nex")
    write_nexus_hmm(
        nexus_path, sequences, network_newick, species_to_samples,
        model=args.model, iterations=args.iterations, threads=args.threads,
        runs=20, outdir=args.outdir
    )
    print(f"[ok] Wrote NEXUS: {nexus_path}")

    # 8) Run PhyloNet
    if not args.write_only:
        cmd = f"{args.phylonet_cmd} {shlex.quote(nexus_path)}"
        run(cmd, check=True, capture=False)
        print("[ok] PhyloNet-HMM finished. Output dir:", args.outdir)

if __name__ == "__main__":
    main()