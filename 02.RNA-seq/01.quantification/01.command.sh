#!/bin/bash
job=$1
## Create an env for sheepGTEx
#conda install -n base -c conda-forge mamba
#mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas

ls 00.mergefq > sample.list
source /storage/public/home/2020060185/anaconda3/envs/snakemake/bin/activate snakemake
#snakemake -s snakefile --configfile config.yaml --dag | dot -Tpdf > dag.pdf
#snakemake -j 4 -s snakefile --configfile config.yaml -p -k --use-conda
snakemake -j 2000 --cluster-config cluster.yaml \
    --cluster "jsub -q {cluster.partition} -M {cluster.mem} -n {cluster.cpus} -R 'span[hosts=1]' -J {cluster.name} -o {cluster.output} -e {cluster.error}" \
    -s snakefile \
    --configfile config.yaml \
    --latency-wait 120 --printshellcmds --keep-going \
    --nolock --rerun-incomplete > sheep.${job}.log 2>&1 
