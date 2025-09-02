#!/bin/bash
## calculate coverage of autosome, chrX and chrY
mkdir -p log indivstat
cat bam.list | while read samfile smid
do
	jsub -q fat -M 30000000 -n 1 -R "span[hosts=1]" -J ${smid}_bam.dp \
		-e log/${smid}.log -o log/${smid}.log \
		"bash 01.indivstat.sh $samfile $smid"
done

