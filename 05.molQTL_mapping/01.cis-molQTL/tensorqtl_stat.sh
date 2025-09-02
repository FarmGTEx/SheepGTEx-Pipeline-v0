#!/bin/bash
tis=$1

# count the number and distance of cis-QTLs
for CHR in chr{1..26}
do
    zcat ${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.${CHR}.txt.gz | awk 'NR!=1{print $1"\t"$2"\t"$3"\t"$7}'
done  | shuf -n 20000 > ${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.shuf20k.txt
