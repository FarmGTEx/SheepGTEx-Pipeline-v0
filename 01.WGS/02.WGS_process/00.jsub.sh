#!/bin/bash
# merge fastq files, input files in the format of: sampleID fastq1  fastq2
mkdir -p Logs/mergefq
for sample in `cut -f1 sheep.mergefq.list | sort -u`
do
        jsub -q normal -M 8000000 -n 1 -R "span[hosts=1]" -J ${sample}.mergefq \
                -e Logs/mergefq/${sample}.%J.log -o Logs/mergefq/${sample}.%J.log \
                "bash 00.mergefq.sh sheep.mergefq.list 00.mergefq $sample"
done