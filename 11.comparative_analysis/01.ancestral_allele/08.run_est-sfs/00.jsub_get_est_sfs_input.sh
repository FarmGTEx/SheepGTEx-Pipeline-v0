#!/bin/bash
for chrom in `cat chr.list`
do
	mkdir chr${chrom}	
	for num in {00..29}
	do
		jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_est_chr${chrom}_${num} -e log/est_chr${chrom}_${num}.%J.log -o log/est_chr${chrom}_${num}.%J.log "bash 00.run_est_sfs.sh ${chrom} ${num}"
	done
done
