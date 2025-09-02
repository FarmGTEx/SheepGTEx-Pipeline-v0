#!/bin/bash
for chrom in chr{1..26} chrX
do
	jsub -q normal -M 8000000 -n 1 -R "span[hosts=1]" -J sheep.beagle.maf.${chrom} \
		-e log/${chrom}.maf.log -o log/${chrom}.maf.log bash 03.maf.sh $chrom
done
