#!/bin/bash
# remove batch effect within each tissue
awk '{print $1"\t"$1"\t"$2}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/sample.list > sample_to_participant.list
mkdir -p tissue
for tis in `cut -f1 tissue.list`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J rmbatch_${tis}_gene \
		-e log/rmbatch_${tis}.%J.log -o log/rmbatch_${tis}.%J.log "bash removeBatchEffect.sh $tis"
done
# merge expression of tissues
jsub -q normal -n 1 -R "span[hosts=1]" -J merge_rmbatch_gene \
	-e log/merge_rmbatch.%J.log -o log/merge_rmbatch.%J.log "bash merge.rmbatch.sh"

# Variance Partition
jsub -q normal -n 3 -R "span[hosts=1]" -J variancePartition_gene_log \
	-e log/variancePartition.log.%J.log -o log/variancePartition.%J.log \
	"Rscript variancePartition.log.R sheep.PCGlnc.gene.merged.tpm.txt info.txt variancePartition.log"
jsub -q normal -n 3 -R "span[hosts=1]" -J variancePartition_rmbatch_gene_log \
	-e log/variancePartition.rmbatch.log.%J.log -o log/variancePartition.rmbatch.log.%J.log \
	"Rscript variancePartition.log.R sheep.PCGlnc.gene.merged.rmbatch.tpm.txt info.txt variancePartition.rmbatch.log"
