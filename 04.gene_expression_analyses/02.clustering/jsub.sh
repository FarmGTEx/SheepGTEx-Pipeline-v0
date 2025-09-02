#!/bin/bash
jsub -q normal -n 4 -R "span[hosts=1]" -J sheep.gene.hclust -e hclust.%J.log -o hclust.%J.log "Rscript hclust.R"
jsub -q normal -n 4 -R "span[hosts=1]" -J sheep.gene.tsne -e tsne.%J.log -o tsne.%J.log "Rscript tsne.R"
