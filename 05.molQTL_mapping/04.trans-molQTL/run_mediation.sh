#!/bin/bash
tis=$1

python mediation.py \
	--info omiga_trans_cis.txt --tissue ${tis} \
	--cis_expr /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/phenotypes/${tis}.expression.bed.gz \
	--trans_expr /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v8.trans/${tis}/phenotypes/${tis}.expression.bed.gz \
	--vcf /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis}.vcf.gz \
	--covar /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/covFile/${tis}.tsv \
	--out /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v8.trans/${tis}/results/omiga/trans_lmm/mediation.csv
