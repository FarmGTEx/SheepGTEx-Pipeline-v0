#!/bin/bash
# prepare metadata (individual, sex and age)
## NOTE: We need assign tissue samples into individuals first based on IBS and GRM matrices, if we are not sure about the individual relationship from public data. (See 05.molQTL_mapping/00.sample_deduplication.ipynb)
awk 'NF==3{print $1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/ind_tis2_sample_melt25.tsv | sed '1d' | sort -u > select.txt # individual list
cat select.txt | sort -u | sed '1iSUBJID' | csvtk join -t -f 'SUBJID;SUBJID' - ../meta_raw.txt > meta.txt # header: SUBJID, SEX, AGE

# prepare expression data
for tis in `cut -f1 ../../tissue25.list | sed '1d'`
do
        mkdir -p ${tis}/phenotypes ${tis}/log
        grep -w $tis /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/ind_tis_sample_melt25.tsv | awk '{print $1"_"$2"\t"$3}' > ${tis}/${tis}.list
        cut -f1 ${tis}/${tis}.list > ${tis}/${tis}.samplelist
done
## prepare expression input of QTL mapping
for tis in `cut -f1 ../../tissue25.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J v3.eqtl_prepare_expression_${tis} -e ${tis}/log/02.${tis}_prepare_pheno.%J.log -o ${tis}/log/02.${tis}_prepare_pheno.%J.log "bash qtl_prepare_expression.sh ${tis}"
done
## merge phenotypes for imputation
mkdir -p log
jsub -q normal -n 4 -R "span[hosts=1]" -J v3.eqtl_expression_merge -e log/merge.%J.log -o log/merge.%J.log "bash merge.sh"
## input expression matrix for HYFA
head -1 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v3.min25_split_impute/observed/express.csv > express.csv # file format: individual, tissue, expression matrix of all genes within tissue (inverse normal transformed TMM)
grep -wf select.txt /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v3.min25_split_impute/observed/express.csv >> express.csv
csvtk cut -f2 express.csv | sed -e '1d' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -nk2 -r > tissue.list
