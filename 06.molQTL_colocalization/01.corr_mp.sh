#!/bin/bash
# prepare data
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/log
	jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_mp_${tis} \
		-e ${tis}/log/01.prepare_mp.${tis}.%J.log -o ${tis}/log/01.prepare_mp.${tis}.%J.log \
		"bash prepare_mp.sh ${tis}"
done

# calculate the correlation of different phenotypes in a gene within a tissue
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J corr_mp_within_tissue_${tis} \
		-e ${tis}/log/01.corr_mp_within.${tis}.%J.log -o ${tis}/log/01.corr_mp_within.${tis}.%J.log \
		"bash corr_mp_within_tissue.sh ${tis}"
done
# combine results
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	echo $tis ; bash corr_mp_within_tissue_combine.sh $tis
done
cat */results/mp_within/mp.corr.txt | sed '1iPairs\tphenotypes\tPearson correlation\ttissue' > mp.corr.txt

# calculate the correlation of different isoform/exon expression in a gene within a tissue
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	for qtl in eeQTL isoQTL
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J corr_mp_within_gene_${tis}_${qtl} \
		-e ${tis}/log/01.corr_mp_within.${qtl}.${tis}.%J.log -o ${tis}/log/01.corr_mp_within.${qtl}.${tis}.%J.log \
		"python corr_mp_within_tissue_isoform.py --infile ${tis}/phenotypes/${qtl}.expression.bed.gz --outfile ${tis}/results/mp_within/${qtl}.within_gene.corr.txt"
	done
done
# combine results
for qtl in eeQTL isoQTL
do
	awk 'NR==1||FNR>1' */results/mp_within/${qtl}.within_gene.corr.txt | sed "s/^/${qtl}\t/g" | sed "1s/^${qtl}/qtl/" > ${qtl}.within_gene.corr.txt
done
