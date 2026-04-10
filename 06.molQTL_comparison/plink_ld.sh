#!/bin/bash
tis=$1
qtl=$2

mkdir -p ${tis}/results/plink
# select lead SNP of eGenes if exist in at least one qtl
csvtk join -t -f 'phenotype_id;phenotype_id' ${tis}/qtls/eQTL.permutation.txt ${tis}/qtls/${qtl}.permutation.txt | csvtk cut -t -f 'phenotype_id,variant_id,is_eGene' | awk '$4=="TRUE"||$5=="TRUE"' | sed "s/^/eQTL\t${qtl}\t/g" > ${tis}/results/plink/eQTL_${qtl}.sig.list
echo -e "pairs\ttissue\tgene\tsnp1\tsnp2\tr2\tD'" > ${tis}/results/plink/eQTL_${qtl}.r2
awk '$4!=$5' ${tis}/results/plink/eQTL_${qtl}.sig.list | while read qtl1 qtl2 gene snp1 snp2 sig1 sig2
do
	plink --sheep --bfile ../01.eQTL/v1.min40_split/${tis}/genotypes/${tis} --ld $snp1 $snp2 | grep 'R-sq' | awk -v t=$tis -v g=$gene -v s1=$snp1 -v s2=$snp2 '{print t"\t"g"\t"s1"\t"s2"\t"$3"\t"$6}' | sed -e "s/^/eQTL_${qtl}\t/g" >> ${tis}/results/plink/eQTL_${qtl}.r2
done
