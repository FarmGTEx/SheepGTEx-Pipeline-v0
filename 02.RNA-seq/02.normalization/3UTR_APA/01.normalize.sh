#!/bin/bash
# merge chromosomes
cut -f1 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/04.MP1/3UTRpolya/all_mapped_reads.txt | sed 's/\//\t/g' | cut -f2 > sample.list
tabix --list-chroms /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/06.recalINFO/chrAuto.filtered.vcf.gz > chrAuto.list
jsub -q normal -n 1 -R "span[hosts=1]" -J merge_apa_quant_res_by_chr -e log/merge_by_chr.%J.log -o log/merge_by_chr.%J.log "bash merge_apa_quant_res_by_chr.sh"

# normalize PDUI in all tissues for clustering
paste sample.list sample.list | csvtk join -t -H -f '1;1' - /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/sample.list > sample_to_participant.list
mkdir -p log
jsub -q normal -n 1 -R "span[hosts=1]" -J normalize_PDUI -e log/norm.%J.log -o log/norm.%J.log "bash normalize.all.sh"

# filter and normalize in each tissue for QTL mapping
for tis in `cut -f1 tissue40.list | sed '1d'`
do
        mkdir -p ${tis}/phenotypes ${tis}/log
        jsub -q normal -n 1 -R "span[hosts=1]" -J ${tis}.3UTRapa \
			-e ${tis}/log/01.normalize.%J.log -o ${tis}/log/01.normalize.%J.log \
			"bash normalize.tis.sh ${tis}"
done