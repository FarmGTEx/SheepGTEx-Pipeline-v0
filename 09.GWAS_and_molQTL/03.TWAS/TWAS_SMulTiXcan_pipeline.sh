#!/bin/bash
trait=$1
SPrediXcan_result=$2
gwas_file=$3
snp_covariance=$4
model_path=$5
result_path=$6
SMulTiXcan=$7

###1. Parallel control (18 threads)
fifoFile="smultixcan_fifo"
rm -f ${fifoFile}
mkfifo ${fifoFile}
exec 9<> ${fifoFile}
rm -f ${fifoFile}

for ((i=0;i<18;i++))
do
    echo "" >&9
done

echo "Start SMulTiXcan jobs for ${trait}..."

###2. Prepare result directory
mkdir -p ${result_path}/${trait}

###3. QTL types
qtl_types="eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL"

###4. Run SMulTiXcan for each QTL type
for type in ${qtl_types}
do
read -u9
{
    /storage/public/home/2020060185/anaconda3/envs/imlabtools/bin/python ${SMulTiXcan} \
        --models_folder ${model_path} \
        --models_name_filter "sheepgtex_(.*)_${type}_models_filtered_signif.db" \
        --models_name_pattern "sheepgtex_(.*)_${type}_models_filtered_signif.db" \
        --snp_covariance ${snp_covariance}/${type}.smultixcan_covariances.txt.gz \
        --metaxcan_folder ${SPrediXcan_result}/${trait}/ \
        --metaxcan_filter "${type}.${trait}.(.*).csv" \
        --metaxcan_file_name_parse_pattern "${type}.${trait}.(.*).csv" \
        --gwas_file ${gwas_file}/${trait}/filtered_STDERR_${trait}1.tbl.txt \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --beta_column b \
        --pvalue_column p \
        --keep_non_rsid \
        --model_db_snp_key varID \
        --cutoff_condition_number 30 \
        --verbosity 7 \
        --throw \
        --output ${result_path}/${trait}/${type}.${trait}.SMultixcan.txt

    echo "" >&9
} &
done

###5. Wait all jobs to finish
wait
exec 9>&-
echo "All SMulTiXcan jobs for ${trait} finished successfully!"
