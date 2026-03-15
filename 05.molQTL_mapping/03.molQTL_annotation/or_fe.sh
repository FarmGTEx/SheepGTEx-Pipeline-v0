#!/bin/bash
tis=$1
qtl=$2

python /storage/public/home/2020060185/script/OR_FE.py \
        ${tis}/${qtl} annotation allsnp.txt ${tis}/${qtl} \
        --maf_ld_file chrAuto.filtered.keep.score.ld
