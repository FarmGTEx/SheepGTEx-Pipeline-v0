#!/bin/bash
interval=$1
dp=$2

chrom=`echo ${interval} | cut -d ":" -f 1`
bcftools view 04.filter/${chrom}.marked.vcf.gz -O z -o split/${interval}.marked.vcf.gz ${interval}
# mask sites according to DP and missing
/storage/public/home/2020060185/anaconda3/envs/pysam/bin/python filter_vcf_depths.py \
	--infile split/${interval}.marked.vcf.gz --depthslist wgs.depthlist${dp} | \
bcftools view -i 'F_MISSING<0.1' -M2 -f PASS --threads 4 | \
bcftools +fill-tags -O z -o split/${interval}.filtered.dp${dp}.vcf.gz --threads 4 -- -t MAF,F_MISSING
bcftools index -f split/${interval}.filtered.dp${dp}.vcf.gz --tbi --threads 4
