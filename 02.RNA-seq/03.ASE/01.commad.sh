#!/bin/bash
ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/08.ASE/*/*_phaser.gene_ae.txt > phaser.gene_ae.list
split -dl 1000 phaser.gene_ae.list # splited into nine files
mkdir -p gene_ae log
for num in x{00..08}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_phaser_pop_$num \
		-e log/prepare_${num}.%J.log -o log/prepare_${num}.%J.log "bash prepare_phaser_pop.sh $num"
done

# Aggregates gene-level haplotypic expression measurement files
jsub -q normal -n 4 -R "span[hosts=1]" -J phaser_pop \
	-e log/phaser_pop.%J.log -o log/phaser_pop.%J.log "bash phaser_pop.sh"
