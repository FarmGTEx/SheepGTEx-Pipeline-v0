#!/bin/bash
job=$1
mkdir -p Logs
# SNP calling 
jsub -q fat -M 32000000 -n 1 -R "span[hosts=1]" -J wgs_snakemake.${job} \
		-e Logs/snakemake.${job}.%J.log -o Logs/snakemake.${job}.%J.log \
	       	"bash 01.command.sh ${job}"
