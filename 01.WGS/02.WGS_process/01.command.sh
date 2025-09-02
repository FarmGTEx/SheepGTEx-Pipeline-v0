#!/bin/bash
job=$1
source /storage/public/home/2020060185/anaconda3/envs/snakemake/bin/activate snakemake
# submit into the cluster
snakemake -j 100 -s snakefile --configfile config.yaml --latency-wait 60 \
    --printshellcmds --keep-going --nolock --rerun-incomplete \
    --cluster-config cluster.yaml \
    --cluster "jsub -q {cluster.partition} -M {cluster.mem} -n {cluster.cpus} -R 'span[hosts=1]' -J {cluster.name} -o {cluster.output} -e {cluster.error}" > snakemake.${job}.log 2>&1
