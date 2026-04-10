#!/bin/bash
#According to https://github.com/wuyangf7/OPERA
# stage 0: data preparation
## prepare besd input for xQTLs
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/opera/input ${tis}/opera/log
	# get genotype data
        jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_genotype_${tis} \
			-e ${tis}/opera/log/03.prepare_gt.${tis}.%J.log -o ${tis}/opera/log/03.prepare_gt.${tis}.%J.log \
			"bash prepare_gt.sh ${tis}"
	for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL
	do
		for chr in chr{1..26}
		do
		jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_BESD_${tis}_${qtl}_${chr} \
			-e ${tis}/opera/log/03.prepare_besd.${qtl}.${chr}.%J.log -o ${tis}/opera/log/03.prepare_besd.${qtl}.${chr}.%J.log \
			"bash prepare_besd.sh ${tis} ${qtl} ${chr} 1"
		done
	done
done
rm -f Embryo/opera/input/3aQTL.* # empty
for tis in `cut -f1 tissue40.list | sed '1d' | grep -v 'Embryo'`
do
	for chr in chr{1..26}
	do
		echo -e "${tis}/opera/input/sQTL.${chr}.besd\n${tis}/opera/input/eeQTL.${chr}.besd\n${tis}/opera/input/isoQTL.${chr}.besd\n${tis}/opera/input/stQTL.${chr}.besd\n${tis}/opera/input/3aQTL.${chr}.besd\n${tis}/opera/input/enQTL.${chr}.besd" | sed 's/.besd$//g' > ${tis}/opera/input/besd.${chr}.list
		echo $tis $chr
	done
done
for tis in Embryo
do
	for chr in chr{1..26}
	do
		echo -e "${tis}/opera/input/sQTL.${chr}.besd\n${tis}/opera/input/eeQTL.${chr}.besd\n${tis}/opera/input/isoQTL.${chr}.besd\n${tis}/opera/input/stQTL.${chr}.besd\n${tis}/opera/input/enQTL.${chr}.besd" | sed 's/.besd$//g' > ${tis}/opera/input/besd.${chr}.list
		echo $tis $chr
	done
done

## split egenes
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/phenotypes/fine_mapping ${tis}/opera/split
	echo $tis
done
## prepare GWAS input for eQTL per tissue per gene and run OPERA
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	for num in `ls ${tis}/opera/split/x?? | sed 's/\//\t/g' | cut -f4`
	do
		jsub -q normal -n 1 -R "span[hosts=1]" -J OPERA_${tis}_${num} \
			-e ${tis}/opera/output/03.opera.${tis}.${num}.%J.log -o ${tis}/opera/output/03.opera.${tis}.${num}.%J.log \
			"bash opera.sh ${tis} ${num} 1"
	done
done
## combine results
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J OPERA_${tis}_combine \
		-e ${tis}/opera/log/03.opera_combine.${tis}.%J.log -o ${tis}/opera/log/03.opera_combine.${tis}.%J.log \
		"bash opera_combine.sh ${tis}"
done
awk 'NR==1 || FNR>1' */opera/output/*stage2.all.prop > opera.stage2.all.prop
awk 'NR==1 || FNR>1' */opera/output/*stage2.own.prop > opera.stage2.own.prop
