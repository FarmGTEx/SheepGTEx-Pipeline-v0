#!/bin/bash
chrom=${1}
python conver_vcf_to_est_sfs_input.py \
--vcf_ref /storage/public/home/2023050465/goatGTEx/04.ancestral_state/01.SheepBased/05.refpanel_vcf/Sheep3125.chr${chrom}.BeaglePhase.rename.AN_AC.vcf.gz \
--vcf_outgroup1 ../06.pairwise_align_vcf_rm_redundancy/chr${chrom}_Sheep_Cattle_rm.vcf.gz \
--vcf_outgroup2 ../06.pairwise_align_vcf_rm_redundancy/chr${chrom}_Sheep_Pig_rm.vcf.gz \
--vcf_outgroup3 ../06.pairwise_align_vcf_rm_redundancy/chr${chrom}_Sheep_Human_rm.vcf.gz \
--outfile chr${chrom}_sfs_input.txt
