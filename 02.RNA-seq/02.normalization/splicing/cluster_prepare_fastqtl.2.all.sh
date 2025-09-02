#!/bin/bash
awk '{print $1"\t"$1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/sample.list > All/All.sample_ind.list

# run filtering in all samples
python3 cluster_prepare_fastqtl.2.py \
	All/All.sample_ind.list \
        gene.bed All All \
        --leafcutter_dir /storage/public/home/2020060185/software/leafcutter-master \
        --input_dir All --output_dir All

