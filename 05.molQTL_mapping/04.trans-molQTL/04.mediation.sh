#!/bin/bash
for tis in `awk 'NR>1&&$2>=200{print $1}' ../tissue40.list`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J mediation_${tis} -e ${tis}/log/04.mediation.${tis}.%J.log -o ${tis}/log/04.mediation.${tis}.%J.log "bash run_mediation.sh ${tis}"
done
cat */results/omiga/trans_lmm/mediation.csv | awk -v FS="," 'NR==1||$1!="tissue"' > mediation.csv
