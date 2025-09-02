#!/bin/bash
tis=$1
source /storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/activate tensorQTL

# 4. QTL mapping for all variant-phenotype pairs and calculate the Inflation Factor (λ) for each gene
python3 tensorqtl_lambda.py \
	--phenotype_bed_file ${tis}/phenotypes/${tis}.expression.bed.gz \
	--covariates_file ${tis}/covFile/${tis}.tsv \
	--plink_prefix_path ${tis}/genotypes/${tis} \
	--outfile ${tis}/results/tensorqtl/trans/${tis}.lambda.txt
