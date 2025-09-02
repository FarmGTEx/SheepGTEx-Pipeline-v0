#!/bin/bash
tis1=$1
tis2=$2

# select eGenes if exist in at least one qtl
cut -f2 combinations/${tis1}_${tis2}/egene.joined.permutation.txt | sed '1d' > combinations/${tis1}_${tis2}/egene.list
echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tphenotype1\tphenotype2\tchrom" > combinations/${tis1}_${tis2}/coloc.pph4
for chr in chr{1..26}
do
	Rscript tissue_coloc_pph4.R \
		/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis1}/results/tensorqtl/nominal/${tis1}.cis_qtl_pairs.${chr}.txt.gz \
		/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis2}/results/tensorqtl/nominal/${tis2}.cis_qtl_pairs.${chr}.txt.gz \
		combinations/${tis1}_${tis2}/egene.list \
		${chr} combinations/${tis1}_${tis2}/coloc.pph4
done
