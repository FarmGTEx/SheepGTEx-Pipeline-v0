#!/bin/bash
# 1. SuSiER
## split egenes
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/phenotypes/fine_mapping ; echo $tis
	zcat ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | awk '$NF=="TRUE"{print $1}' | split -dl 200 ; mv x?? ${tis}/phenotypes/fine_mapping
done
## run SuSiER
for tis in `cut -f1 ../tissue40.list | sed '1d'` 
do
        for num in `ls ${tis}/phenotypes/fine_mapping/x?? | sed 's/\//\t/g' | cut -f4`
        do
			jsub -q normal -n 1 -R "span[hosts=1]" -J eQTL_susier_${tis}_${num} \
				-e ${tis}/log/07.susier.${tis}.${num}.%J.log -o ${tis}/log/07.susier.${tis}.${num}.%J.log \
				"bash susier.sh ${tis} ${num}"
        done
done
## rerun SuSiER
ls */phenotypes/fine_mapping/x?? | sed 's/\//\t/g' | cut -f1,4 > fine_map.list
cat fine_map.list | while read tis num ; do ls ${tis}/log/07.susier.${tis}.${num}.*.log ; done > fine_map.log
ls -d */results/susier/x?? | sed 's/\//\t/g' | cut -f1,4 | while read tis num
do
	rm -f ${tis}/log/07.susier.${tis}.${num}.*.log
	jsub -q normal -n 2 -R "span[hosts=1]" -J eQTL_susier_${tis}_${num} \
		-e ${tis}/log/07.susier.${tis}.${num}.%J.log -o ${tis}/log/07.susier.${tis}.${num}.%J.log \
		"bash susier.sh ${tis} ${num}"
done
## combine results
for tis in `cut -f1 ../tissue40.list | sed '1d'` 
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J eQTL_susier_combine_${tis} \
		-e ${tis}/log/07.susier_combine.${tis}.%J.log -o ${tis}/log/07.susier_combine.${tis}.%J.log \
		"bash susier_combine.sh ${tis}"
done
zcat */results/susier/*.susier.credible.gz | awk 'NR==1||$1!="tissue"' | gzip -c > susier.credible.gz
zcat */results/susier/*.susier.credible.sig.gz | awk 'NR==1||$1!="tissue"' | gzip -c > susier.credible.sig.gz
zcat */results/susier/*.susier.credible.sig.lead.gz | awk 'NR==1||$1!="tissue"' | gzip -c > susier.credible.sig.lead.gz
cat */results/susier/*.crediblenum.txt > susier.crediblenum.txt
cat */results/susier/*.causalnum.txt > susier.causalnum.txt
## compare credible sets with independent QTLs
zcat tensorqtl_independent_qtl.txt.gz | csvtk cut -t -f 'tissue,phenotype_id,rank' | sed '1d' | bedtools groupby -g 1,2 -c 3 -o count -i - | csvtk join -t -H -f '1,2;1,2' susier.crediblenum.txt - | awk '$3>0' | sed '1itissue\tphenotype_id\tcs_count\tnum_rank' > susier_independent.txt
Rscript /storage/public/home/2020060185/script/susie_independent.R susier_independent.txt eQTL.susier_independent.pdf
zcat susier.credible.sig.lead.gz | awk -v OFS="\t" '{print $3,$6,$7,$8,$9,$8"/"$9,"SNP","","\t1"}' | sed -e '1d' -e 's/_/:/g' -e 's/chr//g' | sort -u > susier.credible.sig.lead.chip.txt
zcat susier.credible.sig.lead.gz | awk -v OFS="\t" '$5>0.2{print $3,$6,$7,$8,$9,$8"/"$9,"SNP","","\t1"}' | sed -e '1d' -e 's/_/:/g' -e 's/chr//g' | sort -u > susier.credible.sig.lead.pip0.2.chip.txt
