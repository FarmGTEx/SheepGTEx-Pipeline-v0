#!/bin/bash

mkdir -p input/strong log

### 1. output all top significant eQTL information from permutation results to a .txt file
# input file: extract phenotype_id and variant_id columns from all tissues
awk '$NF=="TRUE"{gsub(/_/,"\t",$8);print $2"\t"$8}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/tensorqtl_permutation.txt | sort -u | sort -V -k2,3 | awk 'BEGIN{print "phenotype_id\tvariant_id\tchr\tpos"}{print $1"\t"$2"_"$3"\t"$2"\t"$3}' > input/strong_pairs.combined_signifpairs.txt
#> output file: input/strong_pairs.combined_signifpairs.txt (colnames: phenotype_id,variant_id,chr,pos)

### 2. extract top pairs from nominal results for each tissue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
    ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.chr*.txt.gz > input/${tis}.nominal_files.txt
    # extract_pairs
    jsub -q normal -n 1 -R "span[hosts=1]" -J extract_pairs_${tis} -e log/01.extract_pairs.${tis}.%J.log -o log/01.extract_pairs.${tis}.%J.log \
        "python3 extract_pairs_tjy.py input/${tis}.nominal_files.txt input/strong_pairs.combined_signifpairs.txt ${tis} -o input/strong"
    #> output file: input/strong/*.extracted_pairs.txt.gz
done

### 3. prepare strong SNP-gene pairs for MashR
ls input/strong/*.extracted_pairs.txt.gz > input/strong_pairs_files.txt
# MashR format file (z-score)
python3 mashr_prepare_input.py input/strong_pairs_files.txt strong_pairs -o input --only_zscore
zcat input/strong_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > strong_pairs.MashR_input.txt.gz
