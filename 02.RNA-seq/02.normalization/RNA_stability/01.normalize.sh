#!/bin/bash
# prepare input
ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/04.MP1/RNA_stability/*/*.consExons.tsv > consExons.pathlist
ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/04.MP1/RNA_stability/*/*.introns.tsv > Introns.pathlist
mkdir -p input
for path in `cat consExons.pathlist Introns.pathlist` ; do ln -s $path input ; done

# get the normalized stability results for all samples
awk '{print $1"\t"$1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/sample.list > sample.list
jsub -q fat -n 1 -R "span[hosts=1]" -J normalize_stability -e norm.%J.log -o norm.%J.log "bash normalize.all.sh"
rm -rf input

# get the normalized stability results for each tissue for QTL mapping
tabix --list-chroms /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/06.recalINFO/chrAuto.filtered.vcf.gz > chrAuto.list
for tis in `cut -f1 tissue40.list | sed '1d'`
do
        mkdir -p ${tis}/phenotypes ${tis}/genotypes ${tis}/covFile ${tis}/results ${tis}/log
        grep -w $tis ind_tis_sample_melt40.tsv | awk '{print $1"_"$2"\t"$3}' > ${tis}/${tis}.list
        cut -f1 ${tis}/${tis}.list > ${tis}/${tis}.samplelist
		## run
		jsub -q normal -n 1 -R "span[hosts=1]" -J ${tis}.normalize_stability \
			-e ${tis}/log/01.normalize.%J.log -o ${tis}/log/01.normalize.%J.log "bash normalize.tis.sh ${tis}"
done

