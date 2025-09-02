#!/bin/bash
tis=$1  ### Liver

#module load julia-1.10.3
export OMIGA_PATH="/storage/public/home/2020060185/software/OmiGA-build-v1.0.1.250113_beta/bin"
export PATH="$OMIGA_PATH:$PATH"
for t in {1..20}
do
omiga --mode trans --verbose --debug --threads 1 --multi-task $t 20 \
      --force-double-precision --covariates ${tis}/covFile/${tis}.tsv \
      --genotype ${tis}/genotypes/${tis} --phenotype ${tis}/phenotypes/${tis}.expression.bed.gz \
      --prefix ${tis} --output-dir ${tis}/results/omiga/trans_lmm
done
