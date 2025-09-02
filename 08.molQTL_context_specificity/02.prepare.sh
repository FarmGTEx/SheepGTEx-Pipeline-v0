#!/bin/bash
group0=All
group1=Europe
group2=Central_and_East_Asia

# 0. prepare input samples (Europe vs. Central_and_East_Asia)
## get tissues with more than 40 samples from Europe vs. Central_and_East_Asia
bash qtl_prepare_tissue.sh
## admixture
for tis in `cat tissue.filtered1.list`
do
        jsub -q normal -n 4 -R "span[hosts=1]" -J v1.ancqtl_admixture_${tis} \
			-e ${tis}/log/02.${tis}_admixture.%J.log -o ${tis}/log/02.${tis}_admixture.%J.log \
			"bash admixture.sh ${tis} 4"
done
## get tissues with suitable samples from ${group1} and ${group2}
echo -e "Tissue\t${group0}\t${group1}\t${group2}" > tissue.filtered2.list
for tis in `cat tissue.filtered1.list`
do
	num1=`cat ${tis}/${tis}.${group0}.list | wc -l`
	num2=`cat ${tis}/${tis}.${group1}.list | wc -l`
	num3=`cat ${tis}/${tis}.${group2}.list | wc -l`
	echo -e "${tis}\t${num1}\t${num2}\t${num3}" | awk '$3>=40&&$4>=40' >> tissue.filtered2.list
done

# 1. prepare expression input of QTL mapping
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.ancqtl_prepare_expression_${tis}_${group} \
		-e ${tis}/log/02.${tis}_${group}_prepare_pheno.%J.log -o ${tis}/log/02.${tis}_${group}_prepare_pheno.%J.log \
		"bash qtl_prepare_expression.sh ${tis} ${group}"
	done
done

# 2. prepare genotype input of QTL mapping
## for all samples
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.ancqtl_prepare_genotypes_${tis}_${group0} \
		-e ${tis}/log/02.${tis}_${group0}_prepare_geno.%J.log -o ${tis}/log/02.${tis}_${group0}_prepare_geno.%J.log \
		"bash qtl_prepare_genotypes.combined.sh ${tis}"
done
## for each group
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group1} ${group2}
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.ancqtl_prepare_genotypes_${tis}_${group} \
		-e ${tis}/log/02.${tis}_${group}_prepare_geno.%J.log -o ${tis}/log/02.${tis}_${group}_prepare_geno.%J.log \
		"bash qtl_prepare_genotypes.sh ${tis} ${group}"
	done
done

# 3. prepare covariates of QTL mapping
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v4.ancqtl_prepare_covariates_${tis}_${group} \
		-e ${tis}/log/02.${tis}_${group}_prepare_cov.%J.log -o ${tis}/log/02.${tis}_${group}_prepare_cov.%J.log \
		"bash qtl_prepare_covariates.sh ${tis} ${group}"
	done
done

# expression and genotype PCA estimated
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	size=`cat ${tis}/genotypes/${tis}.${group}.fam | wc -l`
	phenoElbow=`head -1 ${tis}/covFile/${tis}.${group}_pheno_K_elbow.tsv | awk '{print NF}'`
	echo -e "${tis}\t${group}\tElbow\t${phenoElbow}\t${size}"
	done
done > k.txt
