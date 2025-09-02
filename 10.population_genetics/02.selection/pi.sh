#!/bin/bash
chr=$1
pop=$2
win=$3
step=$4

mkdir -p ${pop}
vcftools --gzvcf vcf/${chr}.vcf.gz --keep ${pop}.id \
	--window-pi $win --window-pi-step $step --out ${pop}/${chr}.${win}_${step}
