#!/bin/bash
cut -f1-6 x00.tpm.txt > x00_18.tpm.txt
for tpm in x{00..18}.tpm.txt ; do csvtk join -j 4 -t -f "1,2,3,4,5,6;1,2,3,4,5,6" -H -O --na 0 x00_18.tpm.txt $tpm > x00_18.tpm.tmp.txt ; mv -f x00_18.tpm.tmp.txt x00_18.tpm.txt ; echo $tpm ; done
sed '1iGeneid' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.list | csvtk join -j 4 -t -f '1;1' x00_18.tpm.txt - | cut -f1,7- > sheep.PCGlnc.gene.merged.tpm.txt
gene_num=`cut -f1 sheep.PCGlnc.gene.merged.tpm.txt | sed '1d' | wc -l`
sample_num=`head -1 sheep.PCGlnc.gene.merged.tpm.txt | awk '{print NF-1}'`
echo -e "#1.2\n${gene_num}\t${sample_num}" > sheep.PCGlnc.gene.merged.tpm.gct
sed 's/Geneid/Name/' sheep.PCGlnc.gene.merged.tpm.txt | awk '{ for (i=1; i<=NF; i++) { if (i==1) { printf $i"\t"$i"\t" } else if (i==NF) {printf $i"\n" } else { printf $i"\t" } } }' | sed 's/Name/Description/2' >> sheep.PCGlnc.gene.merged.tpm.gct
