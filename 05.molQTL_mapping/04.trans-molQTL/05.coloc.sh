#!/bin/bash
# mediated cis-eQTL colocalization
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J coloc_${tis} -e ${tis}/log/05.coloc.${tis}.%J.log -o ${tis}/log/05.coloc.${tis}.%J.log "bash coloc_pph4.sh ${tis}"
done

# combine results
awk 'NR==1||FNR>1' */results/coloc/cis_mediat.trans.pph4 > cis_mediat.trans.pph4
