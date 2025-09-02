#!/bin/bash
tis=$1

# get and remove batch effects of pdui within each tissue
head -1 sheep.PCGlnc.gene.merged.tpm.txt | cut -f7- | sed 's/\t/\n/g' | csvtk join -t -H -f '1;1'  - sample_to_participant.list | grep -w "$tis" > tissue/${tis}.list
/storage/public/home/2020060185/anaconda3/envs/PCAForQTL/bin/Rscript removeBatchEffect.R \
    sheep.PCGlnc.gene.merged.tpm.txt tissue/${tis}.list tissue/${tis}.rmbatch.txt
