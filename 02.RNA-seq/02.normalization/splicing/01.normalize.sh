#!/bin/bash
# pipeline referred to: https://github.com/FarmOmics/ChickenGTEx_pilot_phase/tree/main/molQTL%20mapping/splicing_phenotype_preparation
# preapare
ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/07.WASP_STAR/*/*_Aligned.sortedByCoord.out.filtered.bam.junc > junc.list
awk -v OFS="\t" '$3=="exon"{print $1,$4,$5,$7,$10,$10}' /storage/public/home/2020060185/genome/sheep/reference/sheep.gtf | sed -e 's/;//g' -e 's/"//g' | csvtk -t uniq -H -f 1,2,3,4,5,6 | csvtk join -t -H -f '5;1' - /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.list | sed '1i chr\tstart\tend\tstrand\tgene_id\tgene_name' > exon.bed #1-based coordinate, the same as gtf file
awk -v OFS="\t" '$3=="gene"{print $1,$4,$5,$7,$10,$10}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.gtf | sed -e 's/;//g' -e 's/"//g' | csvtk join -t -H -f '5;1' - /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.list > gene.bed #1-based coordinate, the same as gtf file
sort -k2 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/sample.list | bedtools groupby -i - -g 2 -c 1 -o count | sort -nk2 -r > tissue.list

# clustering in all samples
mkdir -p All log
jsub -q fat -n 1 -R "span[hosts=1]" -J leafcutter.cluster_prepare_fastqtl1 \
	-e log/cluster_prepare_fastqtl1.%J.log -o log/cluster_prepare_fastqtl1.%J.log \
	"bash cluster_prepare_fastqtl.1.sh"

# filtering and normalizing in all samples for clustering
jsub -q fat -n 1 -R "span[hosts=1]" -J all.leafcutter.cluster_prepare_fastqtl2 \
	-e log/cluster_prepare_fastqtl2.all.%J.log -o log/cluster_prepare_fastqtl2.all.%J.log \
	"bash cluster_prepare_fastqtl.2.all.sh"

# filter and normalizing in each tissue for QTL mapping
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	jsub -q fat -n 1 -R "span[hosts=1]" -J ${tis}.cluster_prepare_fastqtl2 \
		-e ${tis}/log/01.cluster_prepare_fastqtl2.%J.log -o ${tis}/log/01.cluster_prepare_fastqtl2.%J.log \
		"bash cluster_prepare_fastqtl.2.tis.sh ${tis}"
done
