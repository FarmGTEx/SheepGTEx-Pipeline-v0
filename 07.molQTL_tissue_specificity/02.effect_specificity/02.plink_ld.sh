#!/bin/bash
tis1=$1
tis2=$2

# select lead SNP of eGenes if exist in at least one qtl
echo -e "gene\tsnp1\tsnp2\tr2\tD'" > combinations/${tis1}_${tis2}/plink.r2

awk '$3!=$9{print $2"\t"$3"\t"$9}' combinations/${tis1}_${tis2}/egene.joined.permutation.txt | while read gene snp1 snp2
do
	plink --sheep --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis1}/genotypes/${tis1} --ld $snp1 $snp2 | grep 'R-sq' | awk -v t=$tis1 -v g=$gene -v s1=$snp1 -v s2=$snp2 '{print t"\t"g"\t"s1"\t"s2"\t"$3"\t"$6}' > combinations/${tis1}_${tis2}/${tis1}.${gene}.r2
	plink --sheep --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis2}/genotypes/${tis2} --ld $snp1 $snp2 | grep 'R-sq' | awk -v t=$tis2 -v g=$gene -v s1=$snp1 -v s2=$snp2 '{print t"\t"g"\t"s1"\t"s2"\t"$3"\t"$6}' > combinations/${tis1}_${tis2}/${tis2}.${gene}.r2
	cat combinations/${tis1}_${tis2}/${tis1}.${gene}.r2 combinations/${tis1}_${tis2}/${tis2}.${gene}.r2 | bedtools groupby -i - -g 2,3,4 -c 5,6 -o collapse,collapse >> combinations/${tis1}_${tis2}/plink.r2
	rm -f combinations/${tis1}_${tis2}/${tis1}.${gene}.r2 combinations/${tis1}_${tis2}/${tis2}.${gene}.r2
done
