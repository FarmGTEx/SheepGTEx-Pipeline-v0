#!/bin/bash
tis=$1
group=$2

rm -f ${tis}/results/susier/${group}/${tis}.crediblenum.txt ${tis}/results/susier/${group}/${tis}.causalnum.txt
for num in `ls ${tis}/phenotypes/${group}/x?? | sed 's/\//\t/g' | cut -f4`
do
	echo $num
	if [ $(wc -l < ${tis}/results/susier/${group}/${num}.susier.gz) -gt 0 ]; then
		# number of credible sets for each gene
		zcat ${tis}/results/susier/${group}/${num}.susier.gz | awk -v FS="\t" '(NF>=13){print $1"\t"$2"\t"$4}' | sed '1d' | sort -uV | bedtools groupby -i - -g 1,2 -c 3 -o count | awk '{print $1"\t"$2"\t"$3-1}' >> ${tis}/results/susier/${group}/${tis}.crediblenum.txt
		# number of causal variants for each credible set
		zcat ${tis}/results/susier/${group}/${num}.susier.gz | awk -v FS="\t" '(NF>=13&&$4!=""){print $1"\t"$2"\t"$4}' | sed '1d' | sort | uniq -c >> ${tis}/results/susier/${group}/${tis}.causalnum.txt
	else
		echo "${tis}/results/susier/${group}/${num}.susier.gz is empty"
	fi
done

# combine susier results
zcat ${tis}/results/susier/${group}/x??.susier.gz | awk 'NR==1||(NF>=13&&$1!="tissue")' | gzip -c > ${tis}/results/susier/${group}/${tis}.susier.gz
# extract 95% credible sets
zcat ${tis}/results/susier/${group}/${tis}.susier.gz | awk -v FS="\t" '$4!=""' | gzip -c > ${tis}/results/susier/${group}/${tis}.susier.credible.gz
# extract eQTLs in 95% credible sets
csvtk join -t -f 'pheno_id,SNP;phenotype_id,variant_id' <(zcat ${tis}/results/susier/${group}/${tis}.susier.credible.gz) <(zcat ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.sig.txt.gz | cut -f1,2,7,10) | gzip -c > ${tis}/results/susier/${group}/${tis}.susier.credible.sig.gz
# extract the lead eQTL in each 95% credible set
zcat ${tis}/results/susier/${group}/${tis}.susier.credible.sig.gz | awk '{if (NR==1){print}else{key=$2" "$4;if ($5>max[key]){max[key]=$5;line[key]=$0}}}END{for (k in line){print line[k]}}' | gzip -c > ${tis}/results/susier/${group}/${tis}.susier.credible.sig.lead.gz
