#!/bin/bash
chrom=$1

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' split/Sheep3125.${chrom}.BeaglePhase.rename.AN_AC.vcf.gz | awk '{if($5/$6>0.5){print $1"\t"$2"\t"$3"\t"$4"\t"1-$5/$6}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5/$6}}' > split/Sheep3125.${chrom}.BeaglePhase.rename.AN_AC.maf
