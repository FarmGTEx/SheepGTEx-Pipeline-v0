#!/bin/bash
# ibs (the same sites)
plink --bfile chrAuto --extract rna_overlapped.pos --sheep --distance 'square' 'ibs' --out chrAuto.overlapped
# grm (the same sites)
gcta-1.94.1 --bfile chrAuto --extract rna_overlapped.pos --make-grm-gz --autosome-num 26 --out chrAuto.overlapped
