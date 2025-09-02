#!/bin/bash
chr=$1
pop=$2
win=$3

vcftools --gzvcf vcf/${chr}.vcf.gz --keep ${pop}.id --TajimaD $win --out ${pop}/${chr}.${win}
