#!/bin/bash
tis=$1

mkdir -p parquet_genotypes/tissue/${tis}
path=/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split
for chr in {1..26}
do
	plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered \
		--sheep --keep ${path}/${tis}/genotypes/${tis}.keep --update-ids ${path}/${tis}/genotypes/${tis}.update --keep-allele-order --chr chr${chr} --make-bed \
		--out parquet_genotypes/tissue/${tis}/genotype.chr${chr}
done
