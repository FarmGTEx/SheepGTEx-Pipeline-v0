#!/bin/bash
tis=$1

# combine results
#zcat ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_*.crossmap.txt.gz | awk 'NR==1||$1!="tissue"' | pigz -c > ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.txt.gz

# trans-eGene detection
Rscript trans_eGene_detection.R \
	${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.txt.gz ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.fdr0.05 0.05

