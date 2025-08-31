#!/bin/bash

python hmm.py \
    --tsi ../../tersect-browser/db-data/mongo-data/gp_data_copy/SGN_aer_hom_snps.tsi \
    --fasta ./SL2.50.fa \
    --region SL2.50ch06:10000000-40000000 \
    --samples \
        S_pim_LYC2798 \
        S_lyc_LA2706 S_lyc_LYC1365 S_lyc_LYC2740 \
    --outdir tomato_introgression_analysis \
    --prefix tomato_chr01_10M \
    --threads 8 \
    --phylonet-cmd "java -jar ../PhyloNet.jar"