#!/bin/bash
tissue=$1

grep -w $tissue ind_tis_sample_melt40.tsv | awk '{print $3"\t"$1"_"$2}' > ${tissue}/${tissue}.sample_ind.list

# run filtering in each tissue
python3 cluster_prepare_fastqtl.2.py \
        ${tissue}/${tissue}.sample_ind.list \
        gene.bed All $tissue \
        --leafcutter_dir /storage/public/home/2020060185/software/leafcutter-master \
        --input_dir All --output_dir ${tissue}/phenotypes

rm -f ${tissue}/phenotypes/${tissue}_perind.counts.filtered.gz.phen_*
rm -f ${tissue}/phenotypes/${tissue}_perind.counts.filtered.gz.qqnorm_*