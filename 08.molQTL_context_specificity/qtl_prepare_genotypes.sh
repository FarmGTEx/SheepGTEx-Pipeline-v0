#!/bin/bash
tis=$1
group=$2

# 2. genotype data
awk '{print $1"\t"$1}' ${tis}/${tis}.${group}.list > ${tis}/genotypes/${tis}.${group}.keep
plink --bfile ${tis}/genotypes/${tis}.All --sheep --keep ${tis}/genotypes/${tis}.${group}.keep \
    --make-bed --keep-allele-order --output-chr chrM --out ${tis}/genotypes/${tis}.${group}
# PCA
plink --bfile ${tis}/genotypes/${tis}.${group} --sheep --pca 2000 'header' 'tabs' --out ${tis}/genotypes/${tis}.${group}
# prepare freq file
plink --bfile ${tis}/genotypes/${tis}.${group} --freq gz --sheep --keep-allele-order --output-chr chrM --out ${tis}/genotypes/${tis}.${group}
