#!/bin/bash
tis=$1

# extract lead eQTL affecting different types of molecular phenotypes
mkdir -p ${tis}/results/pleiotropic/old
mv ${tis}/results/pleiotropic/*.txt ${tis}/results/pleiotropic/*.csv ${tis}/results/pleiotropic/old
for num in {0..3} ; do awk -v t=$tis -v n=$num '$1==t&&$5==n{print $3}' pleiotropic_eqtl.txt > ${tis}/results/pleiotropic/${num}.txt ; done
wc -l ${tis}/results/pleiotropic/*.txt | awk '$1<1{print $2}' | xargs rm

# functional enrichment
python ~/script/OR_FE.py \
	/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/qtl_coloc/${tis}/results/pleiotropic \
	/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/annotation \
	/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/allsnp.txt \
	/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/qtl_coloc/${tis}/results/pleiotropic
