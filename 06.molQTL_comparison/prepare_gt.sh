#!/bin/bash
tis=$1

# prepare genotype
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis} \
    --sheep --make-bed --keep-allele-order --out ${tis}/opera/input/${tis}
