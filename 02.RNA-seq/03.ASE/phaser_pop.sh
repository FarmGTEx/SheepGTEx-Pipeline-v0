#!/bin/bash
#Aggregates gene-level haplotypic expression measurement files across samples to produce a single haplotypic expression matrix, where each row is a gene and each column is a sample
source /storage/public/home/2020060185/anaconda3/envs/phaser/bin/activate phaser
python /storage/public/home/2020060185/software/phaser/phaser_pop/phaser_expr_matrix.py \
    --gene_ae_dir gene_ae \
    --features /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.bed \
    --t 4 --o phaser_expr_matrix
