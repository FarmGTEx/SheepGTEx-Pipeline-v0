#!/bin/bash
mkdir -p log split
for interval in `cat /storage/public/home/2020060185/genome/sheep/reference/sheep.chr.5m.intervals`
do
	for dp in 10
	do 
		awk -v dp=$dp '{print $1"\t"dp"\t"60}' sample.list > wgs.depthlist${dp}
		jsub -q normal -M 32000000 -n 4 -R "span[hosts=1]" -J 02.filter.${interval}.dp${dp} \
			-e log/02.filter.${interval}.dp${dp}.%J.log -o log/02.filter.${interval}.dp${dp}.%J.log \
			"bash 02.filter.sh $interval $dp"
	done
done
