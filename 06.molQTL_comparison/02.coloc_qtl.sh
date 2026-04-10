q1#!/bin/bash
# prepare data
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	echo $tis  ; bash prepare_qtl.sh ${tis}
done

# 0.classify molGenes
for tis in `cut -f1 tissue40.list | sed '1d'`
do
       jsub -q normal -n 1 -R "span[hosts=1]" -J gene_${tis} \
	   	-e ${tis}/log/02.gene.${tis}.%J.log -o ${tis}/log/02.gene.${tis}.%J.log "bash egene_classify.sh ${tis}"
done
awk 'NR==1||FNR>1' */results/isegene.txt > isegene.txt

# 1. calculate ld between eQTL and others
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/log
	for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL
	do
	jsub -q normal -n 1 -R "span[hosts=1]" -J ld_${tis}_${qtl} \
		-e ${tis}/log/02.ld.${tis}.${qtl}.%J.log -o ${tis}/log/02.ld.${tis}.${qtl}.%J.log "bash plink_ld.sh ${tis} ${qtl}"
	done
done
## combine results
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	echo $tis
	cat ${tis}/results/plink/eQTL_*.sig.list | sed "s/^/$tis\t/g" > ${tis}/results/plink/eQTL.sig.list
	cat ${tis}/results/plink/eQTL_*QTL.r2 | grep -v '^pairs' > ${tis}/results/plink/eQTL.r2
done
cat */results/plink/eQTL.sig.list | sed '1itissue\tqtl1\tqtl2\tgene\tsnp1\tsnp2\tis_eGene1\tis_eGene2' > eQTL.sig.list
cat */results/plink/eQTL.r2 | sed "1ipairs\ttissue\tgene\tsnp1\tsnp2\tr2\tD'" > eQTL.r2

# 2. calculate PPH4 between eQTL and others
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/log ${tis}/results/coloc
	for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL
	do
	jsub -q fat -n 1 -R "span[hosts=1]" -J pph4_${tis}_${qtl} \
		-e ${tis}/log/02.pph4.${tis}.${qtl}.%J.log -o ${tis}/log/02.pph4.${tis}.${qtl}.%J.log "bash coloc_pph4.sh ${tis} ${qtl}"
	done
done
## combine results
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	bash coloc_pph4.combine.sh ${tis}
done
awk 'NR==1||FNR>1' */results/coloc/eQTL.addcol.pph4 > eQTL.addcol.pph4
awk 'NR==1||FNR>1' */results/coloc/eQTL.addcol.sig.pph4 > eQTL.addcol.sig.pph4
## functional enrichment of the pleiotropic eQTLs
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/results/pleiotropic
	jsub -q normal -n 1 -R "span[hosts=1]" -J or_fe_${tis} \
		-e ${tis}/log/02.or_fe.%J.log -o ${tis}/log/02.or_fe.%J.log "bash or_fe.sh ${tis}"
done
