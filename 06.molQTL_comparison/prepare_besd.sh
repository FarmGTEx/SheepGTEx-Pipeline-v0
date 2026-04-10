#!/bin/bash
tis=$1
qtl=$2
chr=$3
threads=$4

# stage 0: data preparation
## Make a BESD file from tensorQTL output
if [ ${qtl} != "enQTL" ]
then
	zcat ${tis}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz | awk -v q=$qtl -v OFS="\t" 'NF==9{print q":"$1,$2,$3,$7,$8}' | sed '1d' | pigz > ${tis}/opera/input/${qtl}.${chr}.txt.gz
else ### PCGlnc only for enQTL
	zcat ${tis}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz | sed -e '1s/^/gene_id\t/' -e 's/:/\t/' | csvtk join -t -f 'gene_id;gene_id' - <(sed '1igene_id' /storage/public/home/2020060185/genome/sheep/reference/gene.list) | awk -v q=$qtl -v OFS="\t" 'NF==10{print q":"$1":"$2,$3,$4,$8,$9}' | sed '1d' | pigz > ${tis}/opera/input/${qtl}.${chr}.txt.gz
	
fi
smr --thread-num ${threads} --eqtl-summary ${tis}/opera/input/${qtl}.${chr}.txt.gz --fastqtl-nominal-format --make-besd --out ${tis}/opera/input/${qtl}.${chr}
echo "BESD done."

## Update coordinates of SNPs and genes, frequency of effect allele
### esi
zcat ${tis}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz | awk 'NF==9{print $2"\t"$4}' | sort -u | grep -v 'variant_id' | csvtk join -t -H -j ${threads} -f '1;2' - /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/chrAuto.filtered.bim | awk -v OFS="\t" '{print $3,$1,$4,$5,$6,$7,$2}' > ${tis}/opera/input/${qtl}.${chr}.esi.update
smr --beqtl-summary ${tis}/opera/input/${qtl}.${chr} --thread-num ${threads} --update-esi ${tis}/opera/input/${qtl}.${chr}.esi.update
echo "esi updated."
### epi
tss=/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/gene_tss.txt
if [ ${qtl} == "eQTL" ] || [ ${qtl} == "stQTL" ]
then
	zcat ${tis}/phenotypes/${qtl}.expression.bed.gz | awk -v c=$chr '$1==c{print $1"\t"$4"\t"$3}' | csvtk join -t -H -j ${threads} -f '2;1' - $tss | awk -v q=$qtl -v OFS="\t" '{print $1,q":"$2,0,$3,$2,$6}' | sed 's/^chr//' > ${tis}/opera/input/${qtl}.${chr}.epi.update
else
	zcat ${tis}/phenotypes/${qtl}.expression.bed.gz | awk -v c=$chr '$1==c{print $1"\t"$4"\t"$3}' | csvtk join -t -H -j ${threads} -f '2;1' - ${tis}/phenotypes/${qtl}.phenotype_groups.txt | csvtk join -t -H -j ${threads} -f '4;1' - $tss | awk -v q=$qtl -v OFS="\t" '{print $1,q":"$2,0,$3,$4,$7}' | sed 's/^chr//' > ${tis}/opera/input/${qtl}.${chr}.epi.update
fi
smr --beqtl-summary ${tis}/opera/input/${qtl}.${chr} --thread-num ${threads} --update-epi ${tis}/opera/input/${qtl}.${chr}.epi.update
echo "epi updated."

rm -f ${tis}/opera/input/${qtl}.${chr}.esi.update ${tis}/opera/input/${qtl}.${chr}.bak.esi ${tis}/opera/input/${qtl}.${chr}.epi.update ${tis}/opera/input/${qtl}.${chr}.bak.epi
