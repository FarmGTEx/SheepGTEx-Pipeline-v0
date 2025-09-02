#!/bin/bash
chrom=$1

# mappability
## Step1. Separate the genome
samtools faidx /storage/public/home/2020060185/genome/sheep/reference/sheep.fna ${chrom} > ${chrom}.fa
## Step2. Building the index
genmap index -F ${chrom}.fa -I ${chrom}_index
## Step3. Computing the mappability
genmap map -K 36 -E 2 -I ${chrom}_index -O ${chrom}_36_2 -t -w -bg
genmap map -K 75 -E 2 -I ${chrom}_index -O ${chrom}_75_2 -t -w -bg

