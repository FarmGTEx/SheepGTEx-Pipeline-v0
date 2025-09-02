#!/bin/bash
tis=$1
threads=$2


# prune genotype files
cut -f2 ${tis}/${tis}.grouplist | awk '{print $1"\t"$1}' > ${tis}/admixture.list 
plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered --sheep --keep ${tis}/admixture.list --maf 0.05 --mac 6 --make-bed --keep-allele-order --out ${tis}/genotypes/${tis}.admixture --threads $threads
# prune
plink --bfile ${tis}/genotypes/${tis}.admixture --sheep --indep-pairwise 50 5 0.2 --out ${tis}/genotypes/${tis}.admixture --threads $threads
plink --bfile ${tis}/genotypes/${tis}.admixture --extract ${tis}/genotypes/${tis}.admixture.prune.in --sheep --make-bed --keep-allele-order --out ${tis}/genotypes/${tis}.admixture.pruned --threads $threads
# run admixture
/storage/public/home/2020060185/software/admixture_linux-1.3.0/admixture -j$threads --cv ${tis}/genotypes/${tis}.admixture.pruned.bed 2
mv ${tis}.admixture.pruned.2.* ${tis}/genotypes
awk '{print $1}' ${tis}/genotypes/${tis}.admixture.pruned.fam | csvtk join -t -H -f '1;2' - ${tis}/${tis}.grouplist | paste - ${tis}/genotypes/${tis}.admixture.pruned.2.Q | sort -k4 | sed 's/ /\t/g' > ${tis}/genotypes/${tis}.admixture.pruned.merged.2.Q
# remove mixed sample (ancestry from the group < 0.5)
rm -f ${tis}/${tis}.*.list
bedtools groupby -i ${tis}/genotypes/${tis}.admixture.pruned.merged.2.Q -g 4 -c 5,6 -o mean | awk '{if ($2>$3){print $1"\t5"}else{print $1"\t6"}}' | while read group kline
do
    awk -v g=$group -v k=$kline '$4==g&&$k>0.5{print $2"\t"$1}' ${tis}/genotypes/${tis}.admixture.pruned.merged.2.Q > ${tis}/${tis}.${group}.list
    cut -f1 ${tis}/${tis}.${group}.list > ${tis}/${tis}.${group}.samplelist
done
cat ${tis}/${tis}.*.list > ${tis}/${tis}.All.list
cut -f1 ${tis}/${tis}.All.list > ${tis}/${tis}.All.samplelist
