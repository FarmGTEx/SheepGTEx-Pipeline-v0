#!/bin/bash
# 1. prepare data
## 1.1. prepare annotation
jsub -q fat -n 1 -R "span[hosts=1]" -J snp.ontology -e snp.ontology.%J.log -o snp.ontology.%J.log \
       "bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/06.recalINFO/chrAuto.filtered.eff.vcf.gz | gzip -c > snp.ontology.gz"
bash prepare_annot.sh

## 1.2. prepare qtl
for tis in `cut -f1 tissue40.list | sed '1d'`
do
       jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_qtl_${tis} \
              -e ${tis}/log/qtl.%J.log -o ${tis}/log/qtl.%J.log "bash prepare_qtl.sh ${tis}"
       jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_indep_qtl_${tis} \
              -e ${tis}/log/indep_qtl.%J.log -o ${tis}/log/indep_qtl.%J.log "bash prepare_indep_qtl.sh ${tis}"
done

# 2. OR/FE enrichment
for tis in `cut -f1 tissue40.list | sed '1d'`
do
       jsub -q normal -n 1 -R "span[hosts=1]" -J or_fe_enrichment_${tis} \
              -e ${tis}/log/or_fe.%J.log -o ${tis}/log/or_fe.%J.log "bash or_fe.sh ${tis}"
done
