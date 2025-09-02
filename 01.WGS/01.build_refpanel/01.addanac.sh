#!/bin/bash
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=14G

chr=$1
mkdir -p split/${chr}
# Add "chr"
ln -s /work/home/fan_fengting/01sheep/GTEx_refpanel/rawVCF/Sheep3125_MAF0.01_Miss0.1.snp_BeaglePhase/Sheep3125.${chr}.BeaglePhase.vcf.gz split/${chr}
zcat split/${chr}/Sheep3125.${chr}.BeaglePhase.vcf.gz | awk '{if ($1!~/#/){print "chr"$0}else{print}}' | bgzip > split/${chr}/Sheep3125.${chr}.BeaglePhase.rename.vcf.gz
bcftools index split/${chr}/Sheep3125.${chr}.BeaglePhase.rename.vcf.gz
# Add AN,AC tag
bcftools +fill-tags split/${chr}/Sheep3125.${chr}.BeaglePhase.rename.vcf.gz -Oz --threads 4 -o split/${chr}/Sheep3125.${chr}.BeaglePhase.rename.AN_AC.vcf.gz -- -t AN,AC
bcftools index split/${chr}/Sheep3125.${chr}.BeaglePhase.rename.AN_AC.vcf.gz
