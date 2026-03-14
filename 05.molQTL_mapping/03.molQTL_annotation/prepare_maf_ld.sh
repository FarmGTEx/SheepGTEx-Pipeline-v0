#!/bin/bash
gcta-1.94.1 --bfile chrAuto.filtered.keep --autosome-num 26 --autosome --ld-score --out chrAuto.filtered.keep --threads 8
awk '{print $1"\t"$4"\t"$8}' chrAuto.filtered.keep.score.ld | sed '1d' > maf_ld.txt
