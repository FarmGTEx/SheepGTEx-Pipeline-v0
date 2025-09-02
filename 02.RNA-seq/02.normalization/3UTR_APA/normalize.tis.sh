#!/bin/bash
tissue=$1

grep -w $tissue ind_tis_sample_melt40.tsv | awk '{print $3"\t"$1"_"$2}' > ${tissue}/${tissue}.sample_ind.list

Rscript normalize_3apa_pdui.R -p Dapars2_res.all_chromosomes.txt \
    -t /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/gene_tss.txt \
    -s ${tissue}/${tissue}.sample_ind.list -o ${tissue}/phenotypes/${tissue}.expression

bgzip -f ${tissue}/phenotypes/${tissue}.expression.bed
tabix -f ${tissue}/phenotypes/${tissue}.expression.bed.gz
