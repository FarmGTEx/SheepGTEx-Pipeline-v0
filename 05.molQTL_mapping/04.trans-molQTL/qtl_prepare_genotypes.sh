#!/bin/bash
tis=$1

awk '{print $1"\t"$4-1"\t"$4"\t"$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis}.bim > ${tis}/genotypes/${tis}.cis.bed
# 2. genotype data
bedtools intersect -a ${tis}/genotypes/${tis}.cis.bed -b snp_mappability1.bed > ${tis}/genotypes/${tis}.filter1.bed
bedtools intersect -a ${tis}/genotypes/${tis}.filter1.bed -b sheep.repeat.bed -wa -v > ${tis}/genotypes/${tis}.filter2.bed
cut -f4 ${tis}/genotypes/${tis}.filter2.bed > ${tis}/genotypes/${tis}.filter2.snplist
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis} \
    --sheep --extract ${tis}/genotypes/${tis}.filter2.snplist --make-bed --keep-allele-order --output-chr chrM --out ${tis}/genotypes/${tis}
