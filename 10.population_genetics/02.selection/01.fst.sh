#!/bin/bash
# run fst between two populations
mkdir -p eur_cea eur_mou cea_mou
ln -s eur_cea cea_eur
for chr in chr{1..26}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J eur_cea_${chr}_fst.50000_10000 \
		-e eur_cea/${chr}.50000_10000.%J.log -o eur_cea/${chr}.50000_10000.%J.log \
		"bash fst.sh $chr eur cea 50000 10000"
done

# combine and filter by top 1%
bash fst.combine.sh eur cea 50000 10000
