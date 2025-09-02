#!/bin/bash
chr=$1

newchr=$(($chr-53)) # for sheep NC number to chr number
zcat All.hg38ToRamb2.bed.gz | grep -w "^NC_0560${chr}" | sort -n -k 2 | bedtools groupby -i - -g 1,2,3 -c 5 -o max | sed "s/^NC_0560${chr}../chr$newchr/g" | pigz > output/chr${newchr}.max.bed.gz

# Gene-level phastCons scores
grep -w "^chr${newchr}" /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.bed > output/chr${newchr}.gene.bed
bedtools map -split -sorted -a output/chr${newchr}.gene.bed -b output/chr${newchr}.max.bed.gz -c 4 -o mean -null "" > output/chr${newchr}.phastcons.gene.txt
bedtools coverage -a output/chr${newchr}.gene.bed -b output/chr${newchr}.max.bed.gz > output/chr${newchr}.phastcons.gene.coverage.txt
echo "Gene-level phastCons scores done."
