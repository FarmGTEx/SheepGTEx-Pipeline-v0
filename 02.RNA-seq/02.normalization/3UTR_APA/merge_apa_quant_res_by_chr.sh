#!/bin/bash
Rscript /storage/public/home/2020060185/software/3aQTL-pipe-1.1/src/merge_apa_quant_res_by_chr.R \
        -d /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/04.MP1/3UTRpolya/all \
        -f all -c chrAuto.list -o Dapars2_res.all_chromosomes.txt \
	-s sample.list # WARNING: sample names should be the same order as input files!
