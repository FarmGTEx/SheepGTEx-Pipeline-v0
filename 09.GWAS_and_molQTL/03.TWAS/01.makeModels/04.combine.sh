#!/bin/bash
qtl=$1
tis=$2

# Build a database and filter the database
Rscript make_filter_dbs.R 26 ${qtl}/${tis}

# Combine covariances
awk 'NR==1||$1!="GENE"' ${qtl}/${tis}/covariances/Model_training_chr*_covariances.txt | pigz > ${qtl}/${tis}/dbs/Model_training_covariances.txt.gz

