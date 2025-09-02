#!/bin/bash
tis=$1
group=$2
source /storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/activate tensorQTL

# 3. Conditionally independent cis-QTL mapping
echo 'Conditionally independent cis-QTL mapping'
python3 -m tensorqtl ${tis}/genotypes/${tis}.${group} ${tis}/phenotypes/${tis}.${group}.expression.bed.gz ${tis}/results/tensorqtl/independent/${tis}.${group} \
    --covariates ${tis}/covFile/${tis}.${group}.tsv \
    --cis_output ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.egenes.txt \
    --mode cis_independent
