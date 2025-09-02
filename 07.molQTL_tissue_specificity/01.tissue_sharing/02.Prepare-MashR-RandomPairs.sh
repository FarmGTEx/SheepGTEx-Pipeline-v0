#!/bin/bash

mkdir -p input/random log

### 1. output random nominal SNP-gene pairs information from nominal results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# random select 100000 SNP-gene pairs in each tissue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J combine_signif_pairs_${tis} -e log/02.combine_signif_pairs.${tis}.%J.log -o log/02.combine_signif_pairs.${tis}.%J.log \
        "zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.chr*.txt.gz | cut -f1,2 | grep -v '^phenotype_id' | shuf -n 100000 > input/random/${tis}.nominal_paires.shuf.txt"
done
sort -u input/random/*.nominal_paires.shuf.txt > input/nominal_paires.shuf.txt
awk -v OFS="\t" '{gsub(/_/, "\t",$2);print}' input/nominal_paires.shuf.txt | sort -V -k2,3 | awk 'BEGIN{print "phenotype_id\tvariant_id\tchr\tpos"}{print $1"\t"$2"_"$3"\t"$2"\t"$3}' > input/nominal_pairs.combined_signifpairs.txt
#> output file: nominal_pairs.combined_signifpairs.txt.gz

### 2. extract pairs from nominal results for each tissue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
    # extract_pairs
    jsub -q normal -n 1 -R "span[hosts=1]" -J extract_pairs_${tis} -e log/02.extract_pairs.${tis}.%J.log -o log/02.extract_pairs.${tis}.%J.log \
        "python3 extract_pairs_tjy.py input/${tis}.nominal_files.txt input/nominal_pairs.combined_signifpairs.txt ${tis} -o input/random"
    #> output file: input/random/*.extracted_pairs.txt.gz
done

### 3. prepare 1M random SNP-gene pairs for MashR
ls input/random/*.extracted_pairs.txt.gz > input/random_pairs_files.txt
# MashR format file (z-score)
jsub -q normal -n 1 -R "span[hosts=1]" -J mashr_prepare_input -e log/02.mashr_prepare_input.%J.log -o log/02.mashr_prepare_input.%J.log \
    "python3 mashr_prepare_input.py input/random_pairs_files.txt random_pairs -o input --only_zscore --dropna --subset 1000000 --seed 9823"
zcat input/random_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > random_pairs.MashR_input.txt.gz
