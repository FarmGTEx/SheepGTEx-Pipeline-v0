#!/bin/bash
/storage/public/home/2020060185/anaconda3/envs/gatk4/bin/java -jar \
	/storage/public/home/2020060185/software/snpEff/snpEff.jar \
	ARS-UI_Ramb_v2.0 -csvStats Sheep3518.maf0.01.eff.csv Sheep3518.maf0.01.vcf | bgzip > Sheep3518.maf0.01.eff.vcf.gz

grep 'EXON\|INTRON\|UTR_3_PRIME\|UTR_5_PRIME\|TRANSCRIPT\|SPLICE_\|INTERGENIC\|DOWNSTREAM\|UPSTREAM' Sheep3518.maf0.01.eff.csv | sed -e 's/ , /\t/g' -e 's/%//g' > Sheep3518.maf0.01.eff.feature
