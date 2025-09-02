#!/bin/bash
chr=$1

# merge and generate multi-tissue genotype files
ls parquet_genotypes/tissue/*/genotype.chr${chr}.bed | sed 's/.bed$//g' > parquet_genotypes/raw/bmerge.chr${chr}.list
plink --merge-list parquet_genotypes/raw/bmerge.chr${chr}.list --sheep --keep-allele-order --make-bed --out parquet_genotypes/raw/genotype.chr${chr}
plink --bfile parquet_genotypes/raw/genotype.chr${chr} --sheep --keep-allele-order --chr chr${chr} --recode A --out parquet_genotypes/raw/genotype.chr${chr}
cut -d" " -f2,7- parquet_genotypes/raw/genotype.chr${chr}.raw | sed -e '1s/_[A-Z] / /g' -e '1s/_[A-Z]$//' -e '1s/phenotype_id/varID/' | csvtk transpose -d" " -T > parquet_genotypes/raw/genotype.chr${chr}.txt

# convert to parquet file format
/storage/public/home/2020060185/anaconda3/envs/imlabtools/bin/python model_training_variable_to_parquet.py \
    -variable_file parquet_genotypes/raw/genotype.chr${chr}.txt \
    -parquet_output parquet_genotypes/genotype.chr${chr}.parquet \
    -parsimony 9

rm -f parquet_genotypes/raw/genotype.chr${chr}.*
