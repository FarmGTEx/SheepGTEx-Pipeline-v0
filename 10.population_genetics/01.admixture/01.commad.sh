## 4.1. run admixture to get the ancestry proportion
bash prepare_admixture.sh
for k in {2..10} ; do jsub -q normal -n 4 -R "span[hosts=1]" -J v1.eqtl_prepare_breed_admixture_${k} \
    -e admixture.${k}.%J.log -o admixture.${k}.%J.log "bash admixture.sh $k 4" ; done
bash merge.sh # merge all the results of k

## 4.2 get dist matrix for the mega input to build a NJ tree
bash plink2meg.sh
