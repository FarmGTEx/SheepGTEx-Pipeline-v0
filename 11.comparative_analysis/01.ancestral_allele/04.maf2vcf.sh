#!/bin/bash
chrom=${1}
species=${2}
jsub -q normal -n 1 -R "span[hosts=1]" -J 04.sheep_maftovcf_chr${chrom}_${species} -e log/04.maftovcf_chr${chrom}_${species}.%J.log -o log/04.maftovcf_chr${chrom}_${species}.%J.log "maffilter param=04.param/02.pairwise_maf2vcf/chr${chrom}_Sheep_${species}.vcf.param"
