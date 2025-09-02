#!/bin/bash
for chrom in chr{1..26} chrX
do
	jsub -q normal -M 32000000 -n 4 -R "span[hosts=1]" -J sheep.chunk.${chrom}_5M_2M \
		-e log/02.chunk.${chrom}.5M_2M.%J.log \
		-o log/02.chunk.${chrom}.5M_2M.%J.log \
	       	"bash 02.chunk.sh $chrom 5 2"
done

for chrom in chrY
do
	jsub -q normal -M 32000000 -n 4 -R "span[hosts=1]" -J sheep.chunk.${chrom}_5M_2M \
                -e log/02.chunk.${chrom}.5M_2M.%J.log \
                -o log/02.chunk.${chrom}.5M_2M.%J.log \
	       	"bash 02.chunk.chrY.sh 5 2"
done
