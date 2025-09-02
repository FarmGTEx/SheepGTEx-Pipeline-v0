#!/bin/bash
Rscript normalize_3apa_pdui.R -p Dapars2_res.all_chromosomes.txt \
    -t /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/gene_tss.txt \
    -s sample_to_participant.list -o Dapars2_res.all_chromosomes
