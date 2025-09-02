#!/bin/bash
# prepare SNPs of autosome and chrX
cut -f1,2 /storage/public/home/2020060185/genome/sheep/refpanel/Sheep3518.maf0.01.vcf | grep -v '#\|chrX\|Y' > auto.sites
cut -f1,2 /storage/public/home/2020060185/genome/sheep/refpanel/Sheep3518.maf0.01.vcf | grep '^chrX' > x.sites
# prepare SNPs of chrY of single-copy genes
cut -f1,2 /storage/public/home/2020060185/genome/sheep/refpanel/Sheep3518.maf0.01.vcf | grep '^chrY' > y.sites
python chrY_gff2gtf.py --ingff exonerate75.final.reverse.gff --outgtf exonerate75.final.reverse.gtf
grep -wf scy.genelist exonerate75.final.reverse.gtf | awk '$3=="transcript"' | awk '{print $1"\t"$4-1"\t"$5}' | sort -V | bedtools merge -i - -d 1500000 > scy.bed
awk '{print $1"\t"$2-1"\t"$2}' y.sites | bedtools intersect -a - -b scy.bed | awk '{print $1"\t"$2}' > scy.sites
