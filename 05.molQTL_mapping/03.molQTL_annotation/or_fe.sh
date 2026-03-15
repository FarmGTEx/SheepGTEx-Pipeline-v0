#!/bin/bash
tis=$1
qtl=$2

python /storage/public/home/2020060185/script/OR_FE.py \
        /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/${tis}/${qtl} \
        /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/annotation \
        /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/allsnp.txt \
        /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/functional_enrichment/${tis}/${qtl} \
        --maf_ld_file chrAuto.filtered.keep.score.ld
