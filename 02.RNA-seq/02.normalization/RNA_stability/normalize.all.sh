#!/bin/bash
# assemble stability bed
python stability_assemble_bed.py --input-dir input --samples sample.list --ref_anno ../gene.gtf --output sheep.unnorm.stability.bed

# normalize stability
python normalize_phenotypes.py --input sheep.unnorm.stability.bed --samples sample.list --output sheep.stability.bed
bgzip sheep.stability.bed
