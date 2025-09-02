#!/bin/bash
mkdir -p log
for chrom in chr{1..26}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J sheep.${chrom}.mappability -e log/${chrom}.mappability.log -o log/${chrom}.mappability.log "bash mappability.sh $chrom"
done
jsub -q fat -n 1 -R "span[hosts=1]" -J sheep.mappability.merge -e log/merge.log -o log/merge.log "bash merge.sh"
