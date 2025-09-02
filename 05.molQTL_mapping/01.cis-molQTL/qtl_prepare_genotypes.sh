#!/bin/bash
tis=$1

# 2. genotype data
awk '{print $1"\t"$1}' ${tis}/${tis}.list > ${tis}/genotypes/${tis}.keep
awk '{print $2"\t"$2"\t"$1"\t"$1}' ${tis}/${tis}.list > ${tis}/genotypes/${tis}.update
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered \
    --sheep --keep ${tis}/genotypes/${tis}.keep --maf 0.05 --mac 6 --update-ids ${tis}/genotypes/${tis}.update \
    --make-bed --keep-allele-order --output-chr chrM --out ${tis}/genotypes/${tis}
plink --bfile ${tis}/genotypes/${tis} --sheep --recode vcf-iid 'bgz' --keep-allele-order --output-chr chrM --out ${tis}/genotypes/${tis}
tabix -f ${tis}/genotypes/${tis}.vcf.gz

# PCA
plink --bfile ${tis}/genotypes/${tis} --sheep --pca 2000 'header' 'tabs' --out ${tis}/genotypes/${tis}
