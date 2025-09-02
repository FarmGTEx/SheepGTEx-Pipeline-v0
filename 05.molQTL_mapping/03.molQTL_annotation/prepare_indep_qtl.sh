#!/bin/bash
tis=$1

mkdir -p ${tis}/rank
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk '$NF==1{print $7}' | sort -u > ${tis}/rank/1.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk '$NF==2{print $7}' | sort -u > ${tis}/rank/2.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk '$NF==3{print $7}' | sort -u > ${tis}/rank/3.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk '$NF>=4{print $7}' | sed '1d' | sort -u > ${tis}/rank/4.txt

mkdir -p ${tis}/rank/indep4
zcat ../functional_enrichment/tensorqtl_independent_qtl.txt.gz | awk -v t=$tis '$1==t&&$20>3&&$NF==1{print $8}' | sort -u > ${tis}/rank/indep4/1.txt
zcat ../functional_enrichment/tensorqtl_independent_qtl.txt.gz | awk -v t=$tis '$1==t&&$20>3&&$NF==2{print $8}' | sort -u > ${tis}/rank/indep4/2.txt
zcat ../functional_enrichment/tensorqtl_independent_qtl.txt.gz | awk -v t=$tis '$1==t&&$20>3&&$NF==3{print $8}' | sort -u > ${tis}/rank/indep4/3.txt
zcat ../functional_enrichment/tensorqtl_independent_qtl.txt.gz | awk -v t=$tis '$1==t&&$20>3&&$NF>=4{print $8}' | sort -u > ${tis}/rank/indep4/4.txt

find ${tis}/rank -type f -empty -delete
