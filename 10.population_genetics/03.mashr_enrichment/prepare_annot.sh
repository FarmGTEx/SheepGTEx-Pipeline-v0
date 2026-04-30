#!/bin/bash

mkdir -p annotation
# prepare selected SNPs in different level
awk 'NR>1{print $1"\t"$2-1"\t"$3"\t"$5"\t"$8}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/06.population/02.selection/01.fst/eur_cea/chrAuto.50000_10000.windowed.weir.fst.filter.z | sort -k4,4g > fst.bed
lines=`wc -l /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/06.population/02.selection/01.fst/eur_cea/chrAuto.50000_10000.windowed.weir.fst.filter.z | awk '{printf "%.0f\n",$1/10}'`
split -l $lines -d fst.bed annotation/fst

for i in {00..09}
do
	awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered.bim | bedtools intersect -a - -b annotation/fst${i} | cut -f4 | sort -u > annotation/fst${i}.txt
done
