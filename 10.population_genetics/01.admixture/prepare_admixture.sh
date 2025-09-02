#!/bin/bash
# prune genotype files
cut -f1 ../../../ind_tis_sample_melt40.tsv | sort -u | awk '{print $1"\t"$1}' > ind.list 
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered \
    --sheep --keep ind.list --maf 0.05 --mac 6 --make-bed --keep-allele-order --out chrAuto.filtered.keep
# prune
plink --bfile chrAuto.filtered.keep --sheep --indep-pairwise 50 5 0.2 --out chrAuto.filtered.keep
plink --bfile chrAuto.filtered.keep --extract chrAuto.filtered.keep.prune.in \
    --sheep --make-bed --keep-allele-order --out chrAuto.filtered.keep.pruned
