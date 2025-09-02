#!/bin/bash
for sample in `cut -f1 mergefq.list | sort -u`
do
	jsub -q normal -M 8000000 -n 1 -R "span[hosts=1]" -J ${sample}.mergefq \
		-e Logs/mergefq/${sample}.%J.log -o Logs/mergefq/${sample}.%J.log \
	       	"bash 00.mergefq.sh mergefq.list 00.mergefq $sample"
done
