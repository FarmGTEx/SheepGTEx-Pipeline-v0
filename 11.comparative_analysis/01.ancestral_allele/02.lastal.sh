#!/bin/bash
## last alignment to get the chain file
## https://gitlab.com/mcfrith/last/-/blob/main/doc/last-cookbook.rst
## http://animal.nwsuaf.edu.cn/code/index.php/RGD/loadByGet?address[]=RGD/Documents/Pipeline.php
species=$1
name=$2
version=$3
#mkdir -p ${species}_${name}_${version}
/storage/public/home/2023050465/bin/last-1454/bin/last-train -P4 --revsym --matsym --gapsym -E0.05 -C2 Sheep_ARS-UI_Ramb_v2.0_GCF_016772045.1.fna.db /storage/public/home/2023050465/goatGTEx/04.ancestral_state/00.data/${species}_${name}_${version}.fna > ${species}_${name}_${version}/sheep.mat

/storage/public/home/2023050465/bin/last-1454/bin/lastal -m50 -E0.05 -C2 -P4 -p ${species}_${name}_${version}/sheep.mat Sheep_ARS-UI_Ramb_v2.0_GCF_016772045.1.fna.db /storage/public/home/2023050465/goatGTEx/04.ancestral_state/00.data/${species}_${name}_${version}.fna > ${species}_${name}_${version}/sheep_many.maf

/storage/public/home/2023050465/bin/last-1454/bin/last-split ${species}_${name}_${version}/sheep_many.maf > ${species}_${name}_${version}/sheep.maf
