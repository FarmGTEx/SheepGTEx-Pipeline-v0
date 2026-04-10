#!/bin/bash
tis=$1

for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL
do
	awk -v q=$qtl 'NR>1{print q"\t"$0}' ${tis}/results/mp_within/eQTL_${qtl}.corr.txt
done > ${tis}/results/mp_within/mp.corr.txt
