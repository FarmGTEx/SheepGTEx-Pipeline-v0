#!/bin/bash
chrom=${1}
num=${2}
/storage/public/apps/software/est-sfs/est-sfs-release-2.04/est-sfs \
config.txt \
/storage/public/home/2023050465/goatGTEx/04.ancestral_state/01.SheepBased/rerun_lastal/07.get_est-sfs_input/split/chr${chrom}_sfs_input_${num}.txt \
seedfile.txt \
chr${chrom}/chr${chrom}_sfs_${num}.txt \
chr${chrom}/chr${chrom}_sfs_pvalue_${num}.txt
