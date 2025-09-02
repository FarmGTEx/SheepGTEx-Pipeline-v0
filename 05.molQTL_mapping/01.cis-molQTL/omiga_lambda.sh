#!/bin/bash
tis=$1  ### Liver
threads=$2

#module load julia-1.10.3
export OMIGA_PATH="/storage/public/home/2020060185/software/OmiGA-build-v1.0.1.250113_beta/bin"
export PATH="$OMIGA_PATH:$PATH"
mkdir -p ${tis}/results/omiga/lambda

omiga --mode gwas --lambda-only --verbose --debug --threads $threads --force-double-precision \
      --covariates ${tis}/covFile/${tis}.tsv --genotype ${tis}/genotypes/${tis} \
      --phenotype ${tis}/phenotypes/${tis}.expression.bed.gz \
      --prefix ${tis} --output-dir ${tis}/results/omiga/lambda
