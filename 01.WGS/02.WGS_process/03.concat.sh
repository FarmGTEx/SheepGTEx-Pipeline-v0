#!/bin/bash
ls -1v split/chr*\:*.filtered.dp10.vcf.gz | grep -v 'chrX\|Y' > chrAuto.vcflist
bcftools concat -f chrAuto.vcflist -O z -o chrAuto.vcf.gz --threads 24
bcftools index chrAuto.vcf.gz
plink --vcf chrAuto.vcf.gz --sheep --keep-allele-order --double-id --make-bed --out chrAuto
mv chrAuto.bim chrAuto.bim.old
awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' chrAuto.bim.old > chrAuto.bim
