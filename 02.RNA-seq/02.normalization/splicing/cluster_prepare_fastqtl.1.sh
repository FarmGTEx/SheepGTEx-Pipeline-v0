#!/bin/bash
# cluster junc for samples in all samples
python3 cluster_prepare_fastqtl.1.py \
        junc.list exon.bed All \
        --min_clu_reads 30 --min_clu_ratio 0.001 --max_intron_len 500000 \
        --leafcutter_dir /storage/public/home/2020060185/software/leafcutter-master \
        --output_dir All

rm -f All/*_Aligned.sortedByCoord.out.filtered.bam.junc.All.sorted.gz
