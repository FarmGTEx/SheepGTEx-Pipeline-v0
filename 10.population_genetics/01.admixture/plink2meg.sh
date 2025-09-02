#!/bin/bash
num=`cat ind.list | wc -l`
plink --bfile chrAuto.filtered.keep \
        --sheep --distance-matrix \
        --out chrAuto.filtered.keep
perl plink.distance.matrix.to.mega.pl \
        chrAuto.filtered.keep.mdist.id \
        chrAuto.filtered.keep.mdist $num \
        chrAuto.filtered.keep
