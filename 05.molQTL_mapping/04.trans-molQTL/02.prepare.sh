#!/bin/bash
# prepare input data
tabix --list-chroms /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/06.recalINFO/chrAuto.filtered.vcf.gz > chrAuto.list
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
       	mkdir -p ${tis}/phenotypes ${tis}/genotypes ${tis}/results ${tis}/log
		grep -w $tis ../ind_tis_sample_melt40.tsv | awk '{print $1"_"$2"\t"$3}' > ${tis}/${tis}.list
		cut -f1 ${tis}/${tis}.list > ${tis}/${tis}.samplelist
        ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/covFile ${tis}/covFile
done

# 0. extract non-repeat regions with mappability=1, and genes with mean mappability>=0.8
awk '$4==1' /storage/public/home/2020060185/genome/sheep/reference/mappability/genmap/chrAuto_75_2.bed > snp_mappability1.bed
awk '{print $5"\t"$6-1"\t"$7}' /storage/public/home/2020060185/genome/sheep/reference/sheep.repeatMasker.out | sed '1,3d' > sheep.repeat.bed 
awk '$2>=0.8' /storage/public/home/2020060185/genome/sheep/reference/mappability/crossmap/gene_mappability/gene_mappability.txt > gene_mappability0.8.txt
cut -f1 gene_mappability0.8.txt | sed '1iGeneid' | csvtk join -t -f 'Geneid;Geneid' - sheep.PCGlnc.gene.merged.tpm.txt > sheep.PCGlnc.gene.merged.tpm.txt
gene_num=`cut -f1 sheep.PCGlnc.gene.merged.tpm.txt | sed '1d' | wc -l`
sample_num=`head -1 sheep.PCGlnc.gene.merged.tpm.txt | awk '{print NF-1}'`
echo -e "#1.2\n${gene_num}\t${sample_num}" > sheep.PCGlnc.gene.merged.tpm.gct
sed 's/Geneid/Name/' sheep.PCGlnc.gene.merged.tpm.txt | awk '{ for (i=1; i<=NF; i++) { if (i==1) { printf $i"\t"$i"\t" } else if (i==NF) {printf $i"\n" } else { printf $i"\t" } } }' | sed 's/Name/Description/2' >> sheep.PCGlnc.gene.merged.tpm.gct
cut -f1 gene_mappability0.8.txt | sed '1iGeneid' | csvtk join -t -f 'Geneid;Geneid' - sheep.PCGlnc.gene.merged.count.txt > sheep.PCGlnc.gene.merged.count.txt
gene_num=`cut -f1 sheep.PCGlnc.gene.merged.count.txt | sed '1d' | wc -l`
sample_num=`head -1 sheep.PCGlnc.gene.merged.count.txt | awk '{print NF-1}'`
echo -e "#1.2\n${gene_num}\t${sample_num}" > sheep.PCGlnc.gene.merged.count.gct
sed 's/Geneid/Name/' sheep.PCGlnc.gene.merged.count.txt | awk '{ for (i=1; i<=NF; i++) { if (i==1) { printf $i"\t"$i"\t" } else if (i==NF) {printf $i"\n" } else { printf $i"\t" } } }' | sed 's/Name/Description/2' >> sheep.PCGlnc.gene.merged.count.gct

# 1. prepare expression input of QTL mapping
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
    jsub -q normal -n 1 -R "span[hosts=1]" -J v8.trans_prepare_expression_${tis} \
		-e ${tis}/log/02.${tis}_prepare_pheno.%J.log -o ${tis}/log/02.${tis}_prepare_pheno.%J.log \
		"bash qtl_prepare_expression.sh ${tis}"
done

# 1. prepare genotype input of QTL mapping
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v8.trans_prepare_genotypes_${tis} \
		-e ${tis}/log/02.${tis}_prepare_geno.%J.log -o ${tis}/log/02.${tis}_prepare_geno.%J.log \
		"bash qtl_prepare_genotypes.sh ${tis}"
done

