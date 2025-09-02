#!/bin/bash
tis=$1
t=$2

# 1. cross-mappable filtering
zcat ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_${t}.txt.gz | head -1 | sed 's/^/tissue\t/' | pigz -c > ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_${t}.crossmap.txt.gz
zcat ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_${t}.txt.gz | sed '1d' | while read gene variant af beta beta_se pval
do
	IFS='_' read -r chr pos <<< "$variant"
	echo -e "chr${chr}\t$((${pos}-1))\t${pos}" > ${tis}/results/omiga/trans_lmm/variant.${t}.bed
	cond1=`awk -v g=$gene '{if($1==g){print $2}else if($2==g){print $1}}' /storage/public/home/2020060185/genome/sheep/reference/mappability/crossmap/cross_mappability.txt | wc -l`
	if [ "$cond1" -eq 0 ]; then
		echo -e "$tis\t$gene\t$variant\t$af\t$beta\t$beta_se\t$pval" | pigz -c >> ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_${t}.crossmap.txt.gz
	else
		cond2=`awk -v g=$gene '{if($1==g){print $2}else if($2==g){print $1}}' /storage/public/home/2020060185/genome/sheep/reference/mappability/crossmap/cross_mappability.txt | csvtk join -t -H -f '1;4' - /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/02.imputation/01.genotype/discovery/main/gene_tss.1m.len.bed | cut -f2-4 | bedtools intersect -a - -b ${tis}/results/omiga/trans_lmm/variant.${t}.bed | wc -l` 
		if [ "$cond2" -eq 0 ]; then
			echo -e "$tis\t$gene\t$variant\t$af\t$beta\t$beta_se\t$pval" | pigz -c >> ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.task_${t}.crossmap.txt.gz
		else
			echo "Cis-SNP of cross-mapped gene of $gene found! Filter the variant $variant"
		fi
	fi
done
rm -f ${tis}/results/omiga/trans_lmm/variant.${t}.bed
