#!/bin/bash
tis=$1
source /storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/activate tensorQTL

# 2.1. nominal cis-QTL mapping: nominal associations for all variant-phenotype pairs
echo 'nominal cis-QTL mapping'
python3 -m tensorqtl ${tis}/genotypes/${tis} ${tis}/phenotypes/${tis}.expression.bed.gz ${tis}/results/tensorqtl/nominal/${tis} \
    --covariates ${tis}/covFile/${tis}.tsv --mode cis_nominal

# 2.2. Convert .parquet to .txt.gz
echo 'Convert .parquet to .txt.gz'
for CHR in chr{1..26}
do
    input="${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.${CHR}.parquet"
    output="${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.${CHR}.txt.gz"
    python3 ~/script/parquet2txt.py $input $output
    if [ -f "${output}" ]; then
        rm ${input}
    fi
done
