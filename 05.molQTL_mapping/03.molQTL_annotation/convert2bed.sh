#!/bin/bash
chr=$1

zcat download/chr${chr}.phastCons100way.wigFix.gz | convert2bed --input=wig --max-mem=7G --output=bed -d - | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$5}' | pigz > bed/chr${chr}.phastCons100way.hg38.bed.gz
