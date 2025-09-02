#!/bin/bash
for chrom in `cat chr.list`
do 
jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_input_chr${chrom} -e log/get_input_chr${chrom}.%J.log -o log/get_input_chr${chrom}.%J.log "bash 00.get_est_sfs_input.sh ${chrom}"
done
