#!/bin/bash

# trans-QTL mapping
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
	jsub -q normal -n 4 -R "span[hosts=1]" -J omiga_trans_${tis} \
		-e ${tis}/log/03.omiga_trans.${tis}.%J.log -o ${tis}/log/03.omiga_trans.${tis}.%J.log \
		"bash omiga_trans.sh ${tis} 4"
done

# cross-mappability filtering
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
	for t in {1..20}
	do
		jsub -q normal -n 1 -R "span[hosts=1]" -J trans_crossmap_${tis}_${t} \
			-e ${tis}/log/03.crossmap.${tis}.${t}.%J.log -o ${tis}/log/03.crossmap.${tis}.${t}.%J.log \
			"bash filter_crossmap.sh ${tis} ${t}"
	done
done

# trans-eGene detection
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J trans_eGene_${tis} \
			-e ${tis}/log/03.trans_eGene.${tis}.%J.log -o ${tis}/log/03.trans_eGene.${tis}.%J.log \
			"bash trans_eGene_detection.sh ${tis}"
done
## eGene numbers
echo -e 'Tissue\tNumber of eGenes\tNumber of tested genes\tProportion of eGenes' > egenes.txt
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
	zcat ${tis}/results/omiga/trans_lmm/${tis}.trans_qtl_pairs.crossmap.fdr0.05.txt.gz | sed '1d' | awk -v tis=$tis '{if ($NF=="TRUE"){t+=1}else{f+=1}}END{print tis"\t"t"\t"t+f"\t"t/(t+f)}'
done >> egenes.txt
csvtk join -t -f 'Tissue;Tissue' egenes.txt ../v1.min40_split/egenes.txt > egenes.with_cis.txt
## combine results of trans-eQTLs
zcat */results/omiga/trans_lmm/*.trans_qtl_pairs.crossmap.fdr0.05.txt.gz | awk 'NR==1||$NF=="TRUE"' > omiga_trans_egenes.txt
zcat */results/omiga/trans_lmm/*.trans_qtl_pairs.crossmap.fdr0.05.sig.txt.gz | awk 'NR==1||$NF=="TRUE"' | pigz -c > omiga_trans.sig.txt.gz
## get hotspot
csvtk join -t -f 'tissue,variant_id;tissue,variant_id' omiga_trans_hotspot.txt omiga_trans_egenes.txt | cut -f4 | sed '1d' | sort -u > trans_hotspot.genelist
