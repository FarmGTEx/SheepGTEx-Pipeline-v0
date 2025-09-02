#!/bin/bash
chr=$1

mkdir -p chunk/chr${chr}
zcat bed/chr${chr}.phastCons100way.hg38.bed.gz | split -l 1000000 -d --suffix-length=3 --additional-suffix=.bed - chunk/chr${chr}/chr${chr}.chunk_
ls chunk/chr${chr}/chr${chr}.chunk_* | xargs pigz
