#!/bin/bash

mkdir -p annotation
# prepare background data
cut -f2 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered.bim > allsnp.txt
# prepare ontologies
for ontology in `cat snp.ontology.list`
do
	echo $ontology
	zgrep -w "$ontology" snp.ontology.gz | awk '{print $1"_"$2}' | sed 's/^chr//g' > annotation/${ontology}.txt
done

# prepare chromatin states
for i in {1..15}
do
        echo E${i}
	awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered.bim | bedtools intersect -a - -b /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/AAGs/E${i}_Gs.bed -wa | cut -f4 > annotation/E${i}.txt
done
cat annotation/E{1..5}.txt annotation/E12.txt | sort -u > promoter.txt
cat annotation/E{6..10}.txt | sort -u > enhancer.txt
cat promoter.txt enhancer.txt | sort -u > functional.txt
