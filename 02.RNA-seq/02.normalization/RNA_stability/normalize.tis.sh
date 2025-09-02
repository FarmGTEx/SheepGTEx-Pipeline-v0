#!/bin/bash
tissue=$1

grep -w $tissue ind_tis_sample_melt40.tsv | awk '{print $3"\t"$1"_"$2}' > ${tissue}/${tissue}.sample_ind.list

# normalize stability
python normalize_phenotypes.py --input sheep.unnorm.stability.bed --samples ${tissue}/${tissue}.sample_ind.list \
    --output ${tissue}/phenotypes/${tissue}.expression.bed --chrs chrAuto.list
bgzip -f ${tissue}/phenotypes/${tissue}.expression.bed
tabix -f ${tissue}/phenotypes/${tissue}.expression.bed.gz
