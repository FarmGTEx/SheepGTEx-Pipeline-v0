#!/bin/bash
window=$1
buffer=$2

mkdir -p bins/${window}M_${buffer}M
# 1. chunk: create imputation chunks, for example in chrY (~30 minutes).
GLIMPSE2_chunk_static \
    --input split/Sheep1416.chrY.BeaglePhase.rename.AN_AC.vcf.gz \
    --region chrY --window-mb $window --buffer-mb $buffer \
    --output bins/${window}M_${buffer}M/chunks.chrY.txt --sequential --threads 4
# 2. split_reference: create a binary reference panel for quick reading time.
while IFS="" read -r LINE || [ -n "$LINE" ];
do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    GLIMPSE2_split_reference_static \
        --reference split/Sheep1416.chrY.BeaglePhase.rename.AN_AC.vcf.gz \
        --input-region ${IRG} --output-region ${ORG} --output bins/${window}M_${buffer}M/split
done < bins/${window}M_${buffer}M/chunks.chrY.txt
