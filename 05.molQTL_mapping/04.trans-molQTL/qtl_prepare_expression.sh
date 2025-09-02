#!/bin/bash
tis=$1

# 1. normalized expression
/storage/public/home/2020060185/anaconda3/envs/tensorQTL/bin/python ~/script/eqtl_prepare_expression.py \
	sheep.PCGlnc.gene.merged.tpm.gct sheep.PCGlnc.gene.merged.count.gct \
	/storage/public/home/2020060185/genome/sheep/reference/sheep.gtf \
	../sample_to_participant.list ${tis}/phenotypes/${tis} --chrs chrAuto.list --sample_ids ${tis}/${tis}.samplelist \
	--tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm_inv
