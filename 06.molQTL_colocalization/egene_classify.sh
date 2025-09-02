#!/bin/bash
tis=$1

csvtk cut -t -f 'phenotype_id,is_eGene' ${tis}/qtls/eQTL.permutation.txt | sed '1s/is_eGene/eQTL/' > ${tis}/results/isegene.txt
for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL 
do
	echo $tis $qtl
	csvtk cut -t -f 'phenotype_id,is_eGene' ${tis}/qtls/${qtl}.permutation.txt | sed "1s/is_eGene/${qtl}/" > ${tis}/results/isegene.tmp1 
	csvtk join -t -O -f 'phenotype_id;phenotype_id' ${tis}/results/isegene.txt ${tis}/results/isegene.tmp1 > ${tis}/results/isegene.tmp2
	mv -f ${tis}/results/isegene.tmp2 ${tis}/results/isegene.txt
done
sed -i -e "s/^/${tis}\t/g" -e "1s/^${tis}/tissue/" ${tis}/results/isegene.txt
rm -f ${tis}/results/isegene.tmp1
