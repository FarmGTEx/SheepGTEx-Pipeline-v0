#!/bin/bash
for chrom in `cat chr.list`
do 
jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_summery_chr${chrom} -e log/00_chr${chrom}.%J.log -o log/00_chr${chrom}.%J.log "bash summary_ancestral_state.sh ${chrom}"
done
