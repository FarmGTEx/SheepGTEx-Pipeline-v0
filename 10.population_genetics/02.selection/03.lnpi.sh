#!/bin/bash
# run pi in each population
mkdir -p eur cea eur_cea cea_eur
for chr in chr{1..26}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J eur_${chr}_pi.50000_10000 -e eur/${chr}.50000_10000.%J.log -o eur/${chr}.50000_10000.%J.log "bash pi.sh $chr eur 50000 10000"
	jsub -q normal -n 1 -R "span[hosts=1]" -J cea_${chr}_pi.50000_10000 -e cea/${chr}.50000_10000.%J.log -o cea/${chr}.50000_10000.%J.log "bash pi.sh $chr cea 50000 10000"
done
# concat
win=50000
step=10000
for chr in chr{1..26} ; do sed '1d' eur/${chr}.${win}_${step}.windowed.pi | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t%.6f\n",$5}' ; done | sed "1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" > eur/chrAuto.${win}_${step}.windowed.pi
for chr in chr{1..26} ; do sed '1d' cea/${chr}.${win}_${step}.windowed.pi | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t%.6f\n",$5}' ; done | sed "1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" > cea/chrAuto.${win}_${step}.windowed.pi

# 2. run lnpi ratio
bash lnpi.sh eur cea 50000 10000
bash lnpi.sh cea eur 50000 10000

