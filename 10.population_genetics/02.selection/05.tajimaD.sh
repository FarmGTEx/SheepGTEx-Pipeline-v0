#!/bin/bash
# run tajimaD in each population
mkdir -p eur cea
for chr in chr{1..26}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J eur_${chr}_tajimaD.10000 -e eur/${chr}.10000.%J.log -o eur/${chr}.10000.%J.log "bash tajimaD.sh $chr eur 10000"
	jsub -q normal -n 1 -R "span[hosts=1]" -J cea_${chr}_tajimaD.10000 -e cea/${chr}.10000.%J.log -o cea/${chr}.10000.%J.log "bash tajimaD.sh $chr cea 10000"
done
