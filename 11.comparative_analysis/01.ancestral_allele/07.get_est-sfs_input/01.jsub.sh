#!/bin/bash
for chrom in `cat chr.list`
do 
jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_split_chr${chrom} -e log/01.split_chr${chrom}.%J.log -o log/01.split_chr${chrom}.%J.log "bash 01.split.sh ${chrom}"
done
