#!/bin/bash
tis=$1
qtl=$2

# select eGenes if exist in at least one qtl
csvtk join -t -f 'phenotype_id;phenotype_id' ${tis}/qtls/eQTL.permutation.txt ${tis}/qtls/${qtl}.permutation.txt | csvtk cut -t -f 'phenotype_id,is_eGene' | awk '$2=="TRUE"||$3=="TRUE"{print $1}' > ${tis}/results/coloc/eQTL_${qtl}.sig.genelist
echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tphenotype1\tphenotype2\tchrom\ttissue" > ${tis}/results/coloc/eQTL_${qtl}.pph4
for chr in chr{1..26}
do
	Rscript qtl_coloc_pph4.R \
		${tis}/qtls/eQTL.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz \
		${tis}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz \
		${tis}/results/coloc/eQTL_${qtl}.sig.genelist \
		${tis} ${chr} ${tis}/results/coloc/eQTL_${qtl}.pph4 \
		${tis}/phenotypes/${qtl}.phenotype_groups.txt
done
