#!/bin/bash
jsub -q normal -n 1 -R "span[hosts=1]" -J kmeans_tpm \
    -e log/04.kmeans_tpm.%J.log -o log/04.kmeans_tpm.%J.log "Rscript kmeans.R"

jsub -q normal -n 1 -R "span[hosts=1]" -J rand_index_2 \
    -e log/04.rand_index_2.%J.log -o log/04.rand_index_2.%J.log "Rscript run_rand_index.R"