#!/bin/bash
# 0. download data
mkdir download ; cd download
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.100way.phastCons/md5sum.txt
for file in `awk '{print $2}' md5sum.txt`
do
	wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.100way.phastCons/$file
done
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToGCF_016772045.2.over.chain.gz
cd ..

# 1. convert phastcons wig to bed
mkdir -p bed log
for chr in {1..22} X Y M
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J wig2bed_${chr} \
		-e log/01.${chr}.%J.log -o log/01.${chr}.%J.log "bash convert2bed.sh $chr"
done

# 2. Split bed files to chunks
for chr in {1..22} X Y M
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J chunk_${chr} \
		-e log/02.${chr}.%J.log -o log/02.${chr}.%J.log "bash split.sh $chr"
done

# 3. lift the coordinates from hg38 to Ramb2
for chunk in `ls chunk/*/*.bed.gz`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J liftover_${chunk} \
		-e ${chunk}.%J.log -o ${chunk}.%J.log "bash liftover.sh $chunk"
done

# 4. merge results
cat chunk/chr*/chr*.chunk_*.bed.gz.Mapped.gz > All.hg38ToRamb2.bed.gz
mkdir -p output
for chr in {54..79}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J merge_${chr} \
		-e log/04.NC_0560${chr}.%J.log -o log/04.NC_0560${chr}.%J.log "bash merge.sh $chr"
done

for chr in chr{1..26} ; do cat output/${chr}.phastcons.gene.txt ; done > chrAuto.phastcons.gene.txt
for chr in chr{1..26} ; do cat output/${chr}.phastcons.gene.coverage.txt ; done > chrAuto.phastcons.gene.coverage.txt
for chr in chr{1..26} ; do cat output/${chr}.phastcons.snp.txt ; done > chrAuto.phastcons.snp.txt
csvtk join -t -H -f '4;1' <(awk 'NF==5' chrAuto.phastcons.gene.txt) <(awk '$NF>0.5{print $4"\t"$NF}' chrAuto.phastcons.gene.coverage.txt) > chrAuto.phastcons.gene.cov0.5.txt
