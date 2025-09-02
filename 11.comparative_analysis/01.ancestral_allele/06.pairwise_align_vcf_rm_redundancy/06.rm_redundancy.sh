#!/bin/bash
chrom=${1}
species=${2}
python intersect_ARG-site_outgroup-vcf.py \
--vcf_arg /storage/public/home/2023050465/goatGTEx/04.ancestral_state/01.SheepBased/05.refpanel_vcf/Sheep3125.chr${chrom}.BeaglePhase.rename.AN_AC.vcf.gz \
--raw_vcf /storage/public/home/2023050465/goatGTEx/04.ancestral_state/01.SheepBased/rerun_lastal/04.pairwise_align_vcf/chr${chrom}_Sheep_${species}.vcf.gz \
--out_vcf chr${chrom}_Sheep_${species}_rm.vcf.gz
