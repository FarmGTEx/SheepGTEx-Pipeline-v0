#!/bin/bash
threads=$1

ls -1v vcf/chr*.gz > vcflist
bcftools concat --file-list vcflist -Oz -o vcf/chrAuto.vcf.gz --threads $threads
bcftools index --threads ${threads} vcf/chrAuto.vcf.gz
