#!/bin/bash
qtl=$1

# generate model db folder
mkdir -p ${qtl}/smulitxcan/dbs
for tis in `cut -f1 tissue40.list | sed '1d'`
do
    ln -s $PWD/${qtl}/${tis}/dbs/sheepgtex_models_filtered_signif.db ${qtl}/smulitxcan/dbs/${tis}_filtered_signif.db
done

# make smulitxcan covariances
/storage/public/home/2020060185/anaconda3/envs/imlabtools/bin/python meta_covariance_for_models.py \
    -parquet_genotype_folder parquet_genotypes \
    -parquet_genotype_pattern "genotype.chr(.*).parquet" \
    -model_db_folder ${qtl}/smulitxcan/dbs \
    -model_db_file_pattern "(.*)_filtered_signif.db" \
    -output smulitxcan/${qtl}.smultixcan_covariances.txt.gz \
    -parsimony 9
