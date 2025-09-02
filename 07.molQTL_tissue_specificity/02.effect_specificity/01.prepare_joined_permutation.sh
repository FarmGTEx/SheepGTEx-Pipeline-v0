#!/bin/bash
tis1=$1
tis2=$2

mkdir -p combinations/${tis1}_${tis2}
awk -v t=$tis1 'NR==1||$1==t' tensorqtl_permutation.txt > combinations/${tis1}_${tis2}/${tis1}.permutation.txt
awk -v t=$tis2 'NR==1||$1==t' tensorqtl_permutation.txt > combinations/${tis1}_${tis2}/${tis2}.permutation.txt
# is eGene in at least one tissue
csvtk join -t -f 'phenotype_id;phenotype_id' combinations/${tis1}_${tis2}/${tis1}.permutation.txt combinations/${tis1}_${tis2}/${tis2}.permutation.txt | awk 'NR==1||($7=="TRUE"||$13=="TRUE")' > combinations/${tis1}_${tis2}/egene.joined.permutation.txt

