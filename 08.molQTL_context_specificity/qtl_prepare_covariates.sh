#!/bin/bash
tis=$1
group=$2

# 3. covariates
# 1) phenotype and genotype PC
echo run PCA
/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/run_PCA.R ${tis}/phenotypes/${tis}.${group}.expression.bed.gz ${tis}/covFile/${tis}.${group} ${tis}/genotypes/${tis}.${group}.eigenvec
# 2) other covariates (convert them into 2-category data and drop colinear covariates)
## other known covariates (i.e. sex and age)
echo other known covariates
sed '1iSample\tOld' ${tis}/${tis}.${group}.list | csvtk join -t -f 'Old;Sample' - ../../../../covariates.txt | cut -f1,3- > ${tis}/covFile/${tis}.${group}.knowncovariates.tsv
## calculate adjusted R2 of all covariates using lm function
echo calculate adjusted R2
Rscript ~/script/adjR2.R ${tis}/covFile/${tis}.${group}_pheno_PCs.tsv ${tis}/covFile/${tis}.${group}_geno_PCs.tsv ${tis}/covFile/${tis}.${group}.knowncovariates.tsv ${tis}/covFile/${tis}.${group}
## drop colinear covariates
echo "drop colinear covariates (use sex and age only)"
csvtk cut -t -f "Sample,Sex,Age" ${tis}/covFile/${tis}.${group}.knowncovariates.tsv > ${tis}/covFile/${tis}.${group}.knowncovariates2.tsv
python ~/script/colinear_covariates.py --outfile ${tis}/covFile/${tis}.${group}.known.tsv ${tis}/covFile/${tis}.${group}.knowncovariates2.tsv

# 3) conbine covariates (compare different strategies)
## elbow phenoPCs + 5 genoPCs + known
size=`cat ${tis}/${tis}.${group}.list | wc -l`
if [ $size -lt 200 ]; then
	echo elbow phenoPCs + 5 genoPCs + known
	cut -f1-6 ${tis}/covFile/${tis}.${group}_geno_PCs.tsv > ${tis}/covFile/${tis}.${group}_geno_5PCs.tsv
	/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/combine_covariates.R ${tis}/covFile/${tis}.${group}_pheno_K_elbow.tsv ${tis}/covFile/${tis}.${group}_geno_5PCs.tsv ${tis}/covFile/${tis}.${group}.known.tsv ${tis}/covFile/${tis}.${group}.tsv
else
	echo elbow phenoPCs + 10 genoPCs + known
	cut -f1-11 ${tis}/covFile/${tis}.${group}_geno_PCs.tsv > ${tis}/covFile/${tis}.${group}_geno_10PCs.tsv
	/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/combine_covariates.R ${tis}/covFile/${tis}.${group}_pheno_K_elbow.tsv ${tis}/covFile/${tis}.${group}_geno_10PCs.tsv ${tis}/covFile/${tis}.${group}.known.tsv ${tis}/covFile/${tis}.${group}.tsv
fi
