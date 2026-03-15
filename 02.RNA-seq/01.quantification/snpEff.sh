#!/bin/bash
# SNP ontology annotation
/storage/public/home/2020060185/anaconda3/envs/gatk4/bin/java -jar \
        /storage/public/home/2020060185/software/snpEff/snpEff.jar \
        ARS-UI_Ramb_v2.0 -csvStats 06.recalINFO/chrAuto.filtered.eff.csv 06.recalINFO/chrAuto.filtered.vcf.gz | \
        bgzip > 06.recalINFO/chrAuto.filtered.eff.vcf.gz

grep 'EXON\|INTRON\|UTR_3_PRIME\|UTR_5_PRIME\|TRANSCRIPT\|SPLICE_\|INTERGENIC\|DOWNSTREAM\|UPSTREAM' 06.recalINFO/chrAuto.filtered.eff.csv | sed -e 's/ , /\t/g' -e 's/%//g' > 06.recalINFO/chrAuto.filtered.eff.feature
