#!/bin/bash
## last alignment to get the chain file
## http://animal.nwsuaf.edu.cn/code/index.php/RGD/loadByGet?address[]=RGD/Documents/Pipeline.php
species=$1
name=$2
version=$3
## 寻找两个基因组一对一的序列，比对到多个位置的序列将被过滤
/storage/public/home/2023050465/bin/last-1454/bin/maf-swap ${species}_${name}_${version}/sheep.maf | /storage/public/home/2023050465/bin/last-1454/bin/last-split | /storage/public/home/2023050465/bin/last-1454/bin/maf-swap | /storage/public/home/2023050465/bin/last-1454/bin/last-split | /storage/public/home/2023050465/bin/last-1454/bin/maf-sort > ${species}_${name}_${version}/sheep.one.maf

## sheep2species chain
/storage/public/home/2023050465/bin/last-1454/bin/maf-convert psl ${species}_${name}_${version}/sheep.one.maf > ${species}_${name}_${version}/sheep.psl
## sheep is reference/target (old) genome and speicies is query (new) genome
/storage/public/home/2023050465/bin/kentUtils-302.1.0/bin/faToTwoBit /storage/public/home/2023050465/goatGTEx/04.ancestral_state/00.data/${species}_${name}_${version}.fna /storage/public/home/2023050465/goatGTEx/04.ancestral_state/00.data/${species}_${name}_${version}.2bit
/storage/public/home/2023050465/bin/kentUtils-302.1.0/bin/axtChain -linearGap=medium -psl ${species}_${name}_${version}/sheep.psl Sheep_ARS-UI_Ramb_v2.0_GCF_016772045.1.2bit /storage/public/home/2023050465/goatGTEx/04.ancestral_state/00.data/${species}_${name}_${version}.2bit ${species}_${name}_${version}/sheep2${species}.chain

## split maf by chromosome
/storage/public/home/2023050465/bin/perl-5.38.0/bin/perl maf.rename.species.S.pl ${species}_${name}_${version}/sheep.one.maf Sheep ${species} ${species}_${name}_${version}/sheep.rename.maf
mkdir -p ${species}_${name}_${version}/split
/storage/public/home/2023050465/bin/kentUtils-302.1.0/bin/mafSplit splits.bed ${species}_${name}_${version}/split/ ${species}_${name}_${version}/sheep.rename.maf -useFullSequenceName -byTarget
for i in ${species}_${name}_${version}/split/NW_*.maf ; do sed '1d' $i ; done | sed '1i ##maf version=1 scoring=last' > ${species}_${name}_${version}/split/chrS.maf
rm -f ${species}_${name}_${version}/split/NW_*.maf
