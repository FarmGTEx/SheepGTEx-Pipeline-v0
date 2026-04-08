#!/bin/bash
tis=$1

mkdir -p ${tis}/results/coloc
# select cis-mediating eGenes and trans-eGenes
awk -F ',' -v t=$tis '$1==t&&$NF=="True"{print $4}' mediation.csv | sort -u > ${tis}/results/coloc/cis_mediat.genelist
zcat ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.fdr0.05.sig.txt.gz | sed '1d' | cut -f2 | sort -u > ${tis}/results/coloc/trans.genelist
echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tcis_gene\ttrans_gene\tchrom\ttissue" > ${tis}/results/coloc/cis_mediat.trans.pph4
for chr in chr{1..26}
do
	~/bin/Rscript trans_coloc_pph4.R \
		/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz \
		${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.txt.gz \
		${tis}/results/coloc/cis_mediat.genelist \
		${tis}/results/coloc/trans.genelist \
		${tis} ${chr} ${tis}/results/coloc/cis_mediat.trans.pph4
done
