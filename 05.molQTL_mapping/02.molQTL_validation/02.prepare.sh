#!/bin/bash
# internal validation
awk '$2>=100' ../tissue40.list > tissue100.list
sed '1d' tissue100.list | while read tis size num
do
       	mkdir -p ${tis}/group1/phenotypes ${tis}/group1/genotypes ${tis}/group1/covFile ${tis}/group1/results ${tis}/group1/log
       	mkdir -p ${tis}/group2/phenotypes ${tis}/group2/genotypes ${tis}/group2/covFile ${tis}/group2/results ${tis}/group2/log
	half=$(($size/2))
	grep -w $tis ../ind_tis_sample_melt40.tsv | awk '{print $1"_"$2"\t"$3}' > ${tis}/${tis}.list
	shuf -n $half ${tis}/${tis}.list > ${tis}/group1/${tis}.list
	cut -f1 ${tis}/group1/${tis}.list > ${tis}/group1/${tis}.samplelist
	grep -vwf ${tis}/group1/${tis}.list ${tis}/${tis}.list > ${tis}/group2/${tis}.list
	cut -f1 ${tis}/group2/${tis}.list > ${tis}/group2/${tis}.samplelist
done

# 1. prepare expression input of QTL mapping
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v2.eqtl_prepare_expression_${tis}_${group} \
		-e ${tis}/${group}/log/02.${tis}_prepare_pheno.%J.log -o ${tis}/${group}/log/02.${tis}_prepare_pheno.%J.log \
		"bash qtl_prepare_expression.sh ${tis} ${group}"
	done
done

# 2. prepare genotype input of QTL mapping
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v2.eqtl_prepare_genotypes_${tis}_${group} \
		-e ${tis}/${group}/log/02.${tis}_prepare_geno.%J.log -o ${tis}/${group}/log/02.${tis}_prepare_geno.%J.log \
		"bash qtl_prepare_genotypes.sh ${tis} ${group}"
	done
done


# 3. prepare covariates of QTL mapping
# 3.1 covariates preparation
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v2.eqtl_prepare_covariates_${tis}_${group} \
		-e ${tis}/${group}/log/02.${tis}_prepare_cov.%J.log -o ${tis}/${group}/log/02.${tis}_prepare_cov.%J.log \
		"bash qtl_prepare_covariates.sh ${tis} ${group}"
	done
done

# expression and genotype PCA estimated
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	size=`cat ${tis}/${group}/genotypes/${tis}.fam | wc -l`
	phenoElbow=`head -1 ${tis}/${group}/covFile/${tis}_pheno_K_elbow.tsv | awk '{print NF}'`
	echo -e "${tis}\t${group}\tExpression\tElbow\t${phenoElbow}\t${size}"
	done
done > k.txt
