#!/bin/bash
tis=$1

# site
zcat ${tis}/opera/output/x??.sites.gz | sed '1iTissue\tGene\tNumber' | pigz -c > ${tis}/opera/output/${tis}.stage2.sites.gz

# prop
zcat ${tis}/opera/output/x??.stage2.prop.gz | awk 'NR==1||$1!="Tissue"' | pigz -c > ${tis}/opera/output/${tis}.stage2.prop.gz

# ppa
for i in {1..6}
do
	echo $i
	zcat ${tis}/opera/output/x??.stage2_${i}_expos_assoc.ppa.gz | awk 'NR==1||$1!="Tissue"' | pigz -c > ${tis}/opera/output/${tis}.stage2_${i}_expos_assoc.ppa.gz
done

# combo
zcat ${tis}/opera/output/x??.stage2_combo.res.gz | pigz -c > ${tis}/opera/output/${tis}.stage2_combo.res.gz

# get the proportion of the eGenes (GWAS/QTL) explained by xQTLs
python /storage/public/home/2020060185/script/opera_prop.py \
       --sitefile ${tis}/opera/output/${tis}.stage2.sites.gz \
       --qtl_order sQTL,eeQTL,isoQTL,stQTL,3aQTL,enQTL \
       --outfile_all ${tis}/opera/output/${tis}.stage2.all.prop \
       --outfile_own ${tis}/opera/output/${tis}.stage2.own.prop \
       ${tis}/opera/output/${tis}.stage2_*_expos_assoc.ppa.gz
