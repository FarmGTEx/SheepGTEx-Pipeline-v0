#!/bin/bash
group0=All
group1=Europe
group2=Central_and_East_Asia

# split egenes
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	mkdir -p ${tis}/phenotypes/${group}
	zcat ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.txt.gz | awk '$NF=="TRUE"{print $1}' | split -dl 200
	mv x?? ${tis}/phenotypes/${group} ; echo $tis $group
	done
done

# 1. SuSiER
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
        for num in `ls ${tis}/phenotypes/${group}/x?? | sed 's/\//\t/g' | cut -f4`
        do
			jsub -q normal -n 1 -R "span[hosts=1]" -J susier_${tis}_${group}_${num} \
				-e ${tis}/log/04.susier.${tis}.${group}.${num}.%J.log -o ${tis}/log/04.susier.${tis}.${group}.${num}.%J.log \
				"bash susier.sh ${tis} ${group} ${num}"
        done
	done
done
## rerun SuSiER
ls -d */results/susier/*/x?? | sed 's/\//\t/g' | cut -f1,4,5 | while read tis group num
do
        rm -f ${tis}/log/04.susier.${tis}.${group}.${num}.*.log
        jsub -q normal -n 2 -R "span[hosts=1]" -J susier_${tis}_${group}_${num} \
			-e ${tis}/log/04.susier.${tis}.${group}.${num}.%J.log -o ${tis}/log/04.susier.${tis}.${group}.${num}.%J.log \
			"bash susier.sh ${tis} ${group} ${num}"
done

## combine results
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J susier_combine_${tis}_${group} \
		-e ${tis}/log/04.susier_combine.${tis}.${group}.%J.log -o ${tis}/log/04.susier_combine.${tis}.${group}.%J.log \
		"bash susier_combine.sh ${tis} ${group}"
	done
done
## combine results
for group in ${group0} ${group1} ${group2} ; do zcat */results/susier/${group}/*.susier.credible.gz | awk -v FS="\t" -v g=$group '{if (NR==1){print "group\t"$0}else if($1!="tissue"){print g"\t"$0}}' ; done | awk 'NR==1||$1!="group"' | gzip -c > susier.credible.gz
for group in ${group0} ${group1} ${group2} ; do zcat */results/susier/${group}/*.susier.credible.sig.gz | awk -v FS="\t" -v g=$group '{if (NR==1){print "group\t"$0}else if($1!="tissue"){print g"\t"$0}}' ; done | awk 'NR==1||$1!="group"' | gzip -c > susier.credible.sig.gz
for group in ${group0} ${group1} ${group2} ; do zcat */results/susier/${group}/*.susier.credible.sig.lead.gz | awk -v FS="\t" -v g=$group '{if (NR==1){print "group\t"$0}else if($1!="tissue"){print g"\t"$0}}' ; done | awk 'NR==1||$1!="group"' | gzip -c > susier.credible.sig.lead.gz
for group in ${group0} ${group1} ${group2} ; do cat */results/susier/${group}/*.crediblenum.txt | sed "s/$/\t${group}/g" ; done > susier.crediblenum.txt
for group in ${group0} ${group1} ${group2} ; do cat */results/susier/${group}/*.causalnum.txt | sed "s/$/\t${group}/g" ; done > susier.causalnum.txt
