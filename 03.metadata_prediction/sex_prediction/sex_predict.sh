#!/bin/bash
num_features=$1
method=$2

## referred to: https://github.com/FarmGTEx/metadata-prediction-v1
python sex_predict.py f_classif $num_features $method sheep.PCGlnc.gene.merged.tpm.txt y_train.txt
