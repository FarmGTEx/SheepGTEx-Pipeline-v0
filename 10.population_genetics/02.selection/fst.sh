#!/bin/bash
chr=$1
pop1=$2
pop2=$3
win=$4
step=$5

mkdir -p ${pop1}_${pop2}
vcftools --gzvcf vcf/${chr}.vcf.gz --weir-fst-pop ${pop1}.id --weir-fst-pop ${pop2}.id \
	--fst-window-size $win --fst-window-step $step --out ${pop1}_${pop2}/${chr}.${win}_${step}
