#!/bin/bash
# merge TPM and count
mkdir -p log
split -dl 500 featureCounts.list # splited into 19 files, column format: path    individual  sample
for num in x{00..18} ; do jsub -q normal -M 8000000 -n 1 -R "span[hosts=1]" -J ${num}.merge.geneTPMCount \
    -e log/${num}.merge.%J.log -o log/${num}.merge.%J.log "bash merge.sh $num" ; done

##merge count
jsub -q normal -M 32000000 -n 4 -R "span[hosts=1]" -J sheep.merge.geneCount \
    -e log/count.merge.%J.log -o log/count.merge.%J.log "bash merge_count.sh"

##merge TPM
jsub -q normal -M 32000000 -n 4 -R "span[hosts=1]" -J sheep.merge.geneTPM \
    -e log/tpm.merge.%J.log -o log/tpm.merge.%J.log "bash merge_tpm.sh"