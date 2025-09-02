#!/bin/sh

###1. Parameters
tissue=$1
trait=$2

###2. Paths
summary_path=$3
tmp_path=$4
snp_path=$5
covariance_path=$6
model_path=$7
result_path=$8
SPrediXcan=$9

###3. QTL types
qtl_types="eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL"

###4. Extract GWAS info
cat ${summary_path}/${trait}/*.txt | \
    awk '{print $2, $4, $5, $9, $7}' > ${tmp_path}/${tissue}.${trait}.new.summary.tmp
sed -n '1p' ${tmp_path}/${tissue}.${trait}.new.summary.tmp > ${tmp_path}/${tissue}.${trait}.colnames.tmp

###5. Match SNPs
awk -F " " 'FNR==NR {a[$2]=$2; next}($1 in a){print $0}' \
    ${snp_path}/${tissue}/genotypes/${tissue}.bim \
    ${tmp_path}/${tissue}.${trait}.new.summary.tmp > ${tmp_path}/${tissue}.${trait}.overlap.tmp
cat ${tmp_path}/${tissue}.${trait}.colnames.tmp ${tmp_path}/${tissue}.${trait}.overlap.tmp \
    > ${tmp_path}/${tissue}.${trait}.overlap.txt
rm ${tmp_path}/${tissue}.${trait}.*.tmp

###6. Run S-PrediXcan
# Create result directory automatically
mkdir -p ${result_path}/${trait}

for type in ${qtl_types}
do
    echo "Run ${type} - ${tissue} - ${trait}"
    /storage/public/home/2020060185/anaconda3/envs/imlabtools/bin/python ${SPrediXcan} \
        --model_db_path ${model_path}/${type}/${tissue}/dbs/sheepgtex_models_filtered_signif.db \
        --covariance ${model_path}/${type}/${tissue}/dbs/Model_training_covariances.txt.gz \
        --gwas_file ${tmp_path}/${tissue}.${trait}.overlap.txt \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --beta_column b \
        --pvalue_column p \
        --keep_non_rsid \
        --overwrite \
        --output_file ${result_path}/${trait}/${type}.${trait}.${tissue}.csv
done

###7. Clean
rm ${tmp_path}/${tissue}.${trait}.overlap.txt
