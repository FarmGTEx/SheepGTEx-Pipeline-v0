#!/bin/bash
chrom=${1}
python summary_ancestral_state.py \
--chrom ${chrom} \
--site_file ../07.get_est-sfs_input/chr${chrom}_sfs_input.txt.site \
--freq_file ../07.get_est-sfs_input/chr${chrom}_sfs_input.txt \
--sfs_file ../08.run_est-sfs/chr${chrom}_sfs_pvalue.txt \
--ref_vcf ../../05.refpanel_vcf/Sheep3125.chr${chrom}.BeaglePhase.rename.AN_AC.vcf.gz \
--prefix_outfile chr${chrom}_ancestral_state
