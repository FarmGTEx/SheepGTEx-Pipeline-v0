#!/bin/bash
##prepare for TMM normalization
awk '{print $1"_"$2"\t"$1"_"$2}' ../ind_tis_sample_melt40.tsv > sample_to_participant.list
tabix --list-chroms /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/06.recalINFO/chrAuto.filtered.vcf.gz > chrAuto.list
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/phenotypes ${tis}/genotypes ${tis}/covFile ${tis}/results ${tis}/log
	grep -w $tis ../ind_tis_sample_melt40.tsv | awk '{print $1"_"$2"\t"$3}' > ${tis}/${tis}.list
	cut -f1 ${tis}/${tis}.list > ${tis}/${tis}.samplelist
done

# 1. prepare expression input of QTL mapping
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.eqtl_prepare_expression_${tis} \
		-e ${tis}/log/02.${tis}_prepare_pheno.%J.log -o ${tis}/log/02.${tis}_prepare_pheno.%J.log \
		"bash qtl_prepare_expression.sh ${tis}"
done

# 2. prepare genotype input of QTL mapping
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.eqtl_prepare_genotypes_${tis} \
	-e ${tis}/log/02.${tis}_prepare_geno.%J.log -o ${tis}/log/02.${tis}_prepare_geno.%J.log \
	"bash qtl_prepare_genotypes.sh ${tis}"
done

# 3. prepare covariates of QTL mapping
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.eqtl_prepare_covariates_${tis} \
		-e ${tis}/log/02.${tis}_prepare_cov.%J.log -o ${tis}/log/02.${tis}_prepare_cov.%J.log \
		"bash qtl_prepare_covariates.sh ${tis}"
done

# expression and genotype PCA estimated
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	size=`cat ${tis}/genotypes/${tis}.fam | wc -l`
	phenoElbow=`head -1 ${tis}/covFile/${tis}_pheno_K_elbow.tsv | awk '{print NF}'`
	echo -e "${tis}\tExpression\tElbow\t${phenoElbow}\t${size}"
done > k.txt
