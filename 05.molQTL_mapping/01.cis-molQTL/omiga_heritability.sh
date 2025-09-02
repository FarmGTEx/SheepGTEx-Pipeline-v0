#!/bin/bash
tis=$1  ### Liver
threads=$2

#module load julia-1.10.3
export OMIGA_PATH="/storage/public/home/2020060185/software/OmiGA-build-v1.0.1.250113_beta/bin"
export PATH="$OMIGA_PATH:$PATH"
# run cis heritability
omiga --mode her_est --h2-model Ac --h2-algo em_aireml --verbose --debug --threads $threads --force-double-precision \
      --covariates ${tis}/covFile/${tis}.tsv \
      --genotype ${tis}/genotypes/${tis} --phenotype ${tis}/phenotypes/${tis}.expression.bed.gz \
      --prefix ${tis}.cis --output-dir ${tis}/results/omiga/heritability
# run trans heritability
omiga --mode her_est --h2-model At --h2-algo em_aireml --verbose --debug --threads $threads --force-double-precision \
      --covariates ${tis}/covFile/${tis}.tsv \
      --genotype ${tis}/genotypes/${tis} --phenotype ${tis}/phenotypes/${tis}.expression.bed.gz \
      --prefix ${tis}.trans --output-dir ${tis}/results/omiga/heritability
