#!/bin/bash
tis=$1
group=$2

# 3. covariates
# 1) phenotype and genotype PC
echo run PCA
/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/run_PCA.R ${tis}/${group}/phenotypes/${tis}.expression.bed.gz ${tis}/${group}/covFile/${tis} ${tis}/${group}/genotypes/${tis}.eigenvec
# 2) other covariates (convert them into 2-category data and drop colinear covariates)
## other known covariates (i.e. sex and age)
echo other known covariates
sed '1iSample\tOld' ${tis}/${group}/${tis}.list | csvtk join -t -f 'Old;Sample' - ../../covariates.txt | cut -f1,3- > ${tis}/${group}/covFile/${tis}.knowncovariates.tsv
## drop colinear covariates
echo "drop colinear covariates (use sex and age only)"
csvtk cut -t -f "Sample,Sex,Age" ${tis}/${group}/covFile/${tis}.knowncovariates.tsv > ${tis}/${group}/covFile/${tis}.knowncovariates2.tsv
python ~/script/colinear_covariates.py --outfile ${tis}/${group}/covFile/${tis}.known.tsv ${tis}/${group}/covFile/${tis}.knowncovariates2.tsv

# 3) conbine covariates (compare different strategies)
## elbow phenoPCs + 5 genoPCs + known
echo elbow phenoPCs + 5 genoPCs + known
cut -f1-6 ${tis}/${group}/covFile/${tis}_geno_PCs.tsv > ${tis}/${group}/covFile/${tis}_geno_5PCs.tsv
/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/combine_covariates.R ${tis}/${group}/covFile/${tis}_pheno_K_elbow.tsv ${tis}/${group}/covFile/${tis}_geno_5PCs.tsv ${tis}/${group}/covFile/${tis}.known.tsv ${tis}/${group}/covFile/${tis}.elbow_genoPC5.tsv
## elbow phenoPCs + 10 genoPCs + known
echo elbow phenoPCs + 10 genoPCs + known
cut -f1-11 ${tis}/${group}/covFile/${tis}_geno_PCs.tsv > ${tis}/${group}/covFile/${tis}_geno_10PCs.tsv
/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript ~/script/combine_covariates.R ${tis}/${group}/covFile/${tis}_pheno_K_elbow.tsv ${tis}/${group}/covFile/${tis}_geno_10PCs.tsv ${tis}/${group}/covFile/${tis}.known.tsv ${tis}/${group}/covFile/${tis}.elbow_genoPC10.tsv
