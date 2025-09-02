#!/bin/bash
ls -1v split/Sheep3125.*.BeaglePhase.rename.AN_AC.vcf.gz | grep -v 'chrX' > chrAuto.vcflist
bcftools concat -f chrAuto.vcflist -Oz -o chrAuto.vcf.gz --threads 64
bcftools index chrAuto.vcf.gz --threads 64
plink --vcf chrAuto.vcf.gz --sheep --double-id --make-bed --out chrAuto --threads 64
mv chrAuto.bim chrAuto.bim.old
awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' chrAuto.bim.old > chrAuto.bim
