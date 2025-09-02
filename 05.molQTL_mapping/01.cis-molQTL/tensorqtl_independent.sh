#!/bin/bash
tis=$1
source /storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/activate tensorQTL

zcat ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | awk '{if (NR==1 || $NF=="TRUE") print}' > ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt
# 3. Conditionally independent cis-QTL mapping
echo 'Conditionally independent cis-QTL mapping'
python3 -m tensorqtl ${tis}/genotypes/${tis} ${tis}/phenotypes/${tis}.expression.bed.gz ${tis}/results/tensorqtl/independent/${tis} \
    --covariates ${tis}/covFile/${tis}.tsv \
    --cis_output ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt \
    --mode cis_independent 
    #--phenotype_groups # NOTE: add this option to group multiple phenotypes of a gene
