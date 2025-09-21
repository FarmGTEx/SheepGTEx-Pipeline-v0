#!/bin/bash
# data preparation.
# NOTE: We need assign tissue samples into individuals first based on IBS and GRM matrices, if we are not sure about the individual relationship from public data. (See 05.molQTL_mapping/00.sample_deduplication.ipynb)
awk 'NF==3{print $1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/ind_tis2_sample_melt25.tsv | sed '1d' | sort -u > select.txt # individual list
cat select.txt | sort -u | sed '1iSUBJID' | csvtk join -t -f 'SUBJID;SUBJID' - ../meta_raw.txt > meta.txt # header: SUBJID, SEX, AGE
head -1 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v3.min25_split_impute/observed/express.csv > express.csv # file format: individual, tissue, expression matrix of all genes within tissue (inverse normal transformed TMM)
grep -wf select.txt /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v3.min25_split_impute/observed/express.csv >> express.csv
csvtk cut -f2 express.csv | sed -e '1d' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -nk2 -r > tissue.list
