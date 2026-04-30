#!/bin/bash
i=$1

python ~/script/OR_FE.py \
	lfsr${i}/tissues annotation /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/allsnp.txt lfsr${i}/tissues \
	--maf_ld_file /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v9.breed_specific/admixture/chrAuto.filtered.keep.score.ld
