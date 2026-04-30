#!/bin/bash
# 1. prepare annotation and qtl
bash prepare_annot.sh
bash prepare_qtl.sh

# 2. OR/FE enrichment
for i in {1..3} {49..51}
do
    jsub -q normal -n 1 -R "span[hosts=1]" -J or_fe_lfsr${i} \
    -e lfsr${i}/log/or_fe.%J.log -o lfsr${i}/log/or_fe.%J.log "bash or_fe.sh ${i}"
done
