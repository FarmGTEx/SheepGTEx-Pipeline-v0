#!/bin/bash
qtl=$1
tis=$2

##Extract the corresponding information columns from the VCF file to generate a whole chromosome SNP annotation file
if [ ${qtl} == "3aQTL" ]
then
	awk -v OFS="\t" '{print $1,$4,$2,$6,$5,$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/genotypes/${tis}.bim | sed 's/^chr//g' > ${qtl}/${tis}/data/snp_annotation.txt
else # the genotypes of other QTLs are the same as eQTL
	awk -v OFS="\t" '{print $1,$4,$2,$6,$5,$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis}.bim | sed 's/^chr//g' > ${qtl}/${tis}/data/snp_annotation.txt
fi
##Prepare corresponding gene expression, covariate and residual files
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/qtl_coloc/${tis}/phenotypes/${qtl}.expression.bed.gz | cut -f4- | csvtk transpose -t | grep -v '^gene_id' > ${qtl}/${tis}/output/transformed_expression.txt
Rscript preprocess_expr.R $qtl $tis
##Processing genotype files by chromosome
for chr in {1..26}
do
	awk -v c=$chr '$1==c' ${qtl}/${tis}/data/snp_annotation.txt | sed '1i\chromosome\tpos\tvarID\tref_vcf\talt_vcf\trsid' > ${qtl}/${tis}/output/snp_annot.chr${chr}.txt
	if [ ${qtl} == "3aQTL" ]
	then
		plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/genotypes/${tis} --sheep --keep-allele-order --chr chr${chr} --recode A --out ${qtl}/${tis}/data/genotype.chr${chr}
	else
		plink --bfile /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/genotypes/${tis} --sheep --keep-allele-order --chr chr${chr} --recode A --out ${qtl}/${tis}/data/genotype.chr${chr} 
	fi
	cut -d" " -f2,7- ${qtl}/${tis}/data/genotype.chr${chr}.raw | csvtk join -d" " -T -f '1;1' <(cut -f1 ${qtl}/${tis}/output/transformed_expression.txt) - | sed -e '1s/_[A-Z]\t/\t/g' -e '1s/_[A-Z]$//' -e '1s/phenotype_id/varID/' | csvtk transpose -t > ${qtl}/${tis}/output/genotype.chr${chr}.txt
	rm -f ${qtl}/${tis}/data/genotype.chr${chr}.raw
done
