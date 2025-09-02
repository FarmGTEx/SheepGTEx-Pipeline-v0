#!/bin/bash
cat indivstat/* | grep -vw 'nan' > bam.out
/storage/public/home/2020060185/anaconda3/envs/pysam/bin/python genderCheck_snpDepth.py popcluster \
	--statfile bam.out --outprefix bam
