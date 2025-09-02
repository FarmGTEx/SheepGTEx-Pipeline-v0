#!/bin/bash
samfile=$1
smid=$2

/storage/public/home/2020060185/anaconda3/envs/pysam/bin/python genderCheck_snpDepth.py indivstat \
	--samfile $samfile --smid $smid \
	--autosites auto.sites --chrxsites x.sites --chrysites scy.sites \
	--outprefix indivstat/$smid 
