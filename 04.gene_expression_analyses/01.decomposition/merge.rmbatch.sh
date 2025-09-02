#!/bin/bash
cut -f1 tissue.list | head -1 | while read tis ; do cp tissue/${tis}.rmbatch.txt sheep.PCGlnc.gene.merged.rmbatch.tpm.txt ; done
for tis in `cut -f1 tissue.list | sed '1d'` ; do echo $tis ; csvtk join -t -C "$" -f '1,2,3,4,5,6;1,2,3,4,5,6' sheep.PCGlnc.gene.merged.rmbatch.tpm.txt tissue/${tis}.rmbatch.txt > tmp.rmbatch.txt ; mv tmp.rmbatch.txt sheep.PCGlnc.gene.merged.rmbatch.tpm.txt ; done
