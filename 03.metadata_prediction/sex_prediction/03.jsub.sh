#!/bin/bash
# prepare expression and training data (true set)
cut -f1,7- /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene/stringtie/sheep.PCGlnc.gene.merged.tpm.txt > sheep.PCGlnc.gene.merged.tpm.txt
# perform sex predicting
fn=20
method=mlp
jsub -q normal -n 5 -R "span[hosts=1]" -J sex_predict_F${fn}_$method \
	-e sex_predict.F${fn}.${method}.%J.log -o sex_predict.F${fn}.${method}.%J.log \
	"bash sex_predict.sh ${fn} ${method}"
