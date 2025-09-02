#!/bin/bash
tis=$1  ### Liver
threads=$2

#module load julia-1.10.3
export OMIGA_PATH="/storage/public/home/2020060185/software/OmiGA-build-v1.0.1.250113_beta/bin"
export PATH="$OMIGA_PATH:$PATH"
omiga --mode cis --verbose --debug --threads $threads --force-double-precision \
      --covariates ${tis}/covFile/${tis}.tsv \
      --genotype ${tis}/genotypes/${tis} --phenotype ${tis}/phenotypes/${tis}.expression.bed.gz \
      --prefix ${tis} --output-dir ${tis}/results/omiga/cis_lmm

# calculate corelation of nominal results between tensorqtl and omiga
python omiga_cis_stat.py --chrlist chrAuto.list --indir ${tis}/results --outfile ${tis}/results/omiga/cis_lmm/omiga.stat.txt
