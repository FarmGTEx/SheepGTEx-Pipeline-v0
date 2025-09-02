#!/bin/bash
for chrom in chr{1..26} ; do cat ${chrom}_36_2.bedgraph ; done > chrAuto_36_2.bed
for chrom in chr{1..26} ; do cat ${chrom}_75_2.bedgraph ; done > chrAuto_75_2.bed
### Prepare the bed file with mappability < 0.5
awk '{if ($4<0.5) print}' chrAuto_75_2.bed | bedtools sort -i - | bedtools merge -i - > chrAuto_75_2_0.5.bed
awk '{if ($4<0.5) print}' chrAuto_36_2.bed | bedtools sort -i - | bedtools merge -i - > chrAuto_36_2_0.5.bed
