#!/bin/bash
tis=$1
group=$2

# 2. genotype data
awk '{print $1"\t"$1}' ${tis}/${group}/${tis}.list > ${tis}/${group}/genotypes/${tis}.keep
awk '{print $2"\t"$2"\t"$1"\t"$1}' ${tis}/${group}/${tis}.list > ${tis}/${group}/genotypes/${tis}.update
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered --sheep --keep ${tis}/${group}/genotypes/${tis}.keep --maf 0.05 --mac 6 --update-ids ${tis}/${group}/genotypes/${tis}.update \
    --make-bed --keep-allele-order --output-chr chrM --out ${tis}/${group}/genotypes/${tis}
# PCA
plink --bfile ${tis}/${group}/genotypes/${tis} --sheep --pca 2000 'header' 'tabs' --out ${tis}/${group}/genotypes/${tis}
