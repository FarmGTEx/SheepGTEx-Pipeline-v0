#!/bin/bash
chr=$1
threads=$2

bcftools view --threads ${threads} -S ../../query.samplelist /storage/public/home/2020060185/genome/sheep/refpanel/chrAuto.vcf.gz ${chr} | \
bcftools view --threads ${threads} -i 'MAF>=0.05 & MAC>=6' -Oz -o vcf/${chr}.vcf.gz
bcftools index --threads ${threads} vcf/${chr}.vcf.gz
