#!/bin/bash
tis=$1
source /storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/activate tensorQTL

# 1.1. permutation cis-QTL mapping: empirical p-values for phenotypes
echo 'permutation cis-QTL mapping'
python3 -m tensorqtl ${tis}/genotypes/${tis} ${tis}/phenotypes/${tis}.expression.bed.gz ${tis}/results/tensorqtl/permutation/${tis} \
    --covariates ${tis}/covFile/${tis}.tsv --mode cis --seed 9527 
    #--phenotype_groups # NOTE: add this option to group multiple phenotypes of a gene

# 1.2. Define eGenes (FDR 0.05)
echo 'Define eGene (FDR 0.05)'
Rscript ~/script/eGene_detection.R ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl.txt.gz \
    ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz 0.05
