#!/bin/bash
for chrom in `cat chr.list`
do 
for species in {'Cattle','Pig','Human'}
do
jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_rm_redundancy_chr${chrom}_${species} -e log/rm_chr${chrom}_${species}.%J.log -o log/rm_chr${chrom}_${species}.%J.log "bash 00.rm_redundancy.sh ${chrom} ${species}"
done
done
