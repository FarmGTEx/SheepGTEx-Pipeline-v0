#!/bin/bash
tis=$1

mkdir -p ${tis}/results/mp_within
# eQTL vs others
#for mp in sQTL eeQTL isoQTL stQTL 3aQTL 
for mp in eeQTL
do
	python corr_mp_within_tissue.py --infile1 ${tis}/phenotypes/eQTL.expression.bed.gz --infile2 ${tis}/phenotypes/${mp}.expression.bed.gz --outfile ${tis}/results/mp_within/eQTL_${mp}.corr.txt
done
