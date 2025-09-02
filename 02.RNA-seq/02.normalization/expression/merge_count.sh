#!/bin/bash
cut -f1-6 x00.count.txt > x00_18.count.txt
for count in x{00..18}.count.txt ; do csvtk join -j 4 -t -f "1,2,3,4,5,6;1,2,3,4,5,6" -H -O --na 0 x00_18.count.txt $count > x00_18.count.tmp.txt ; mv -f x00_18.count.tmp.txt x00_18.count.txt ; echo $count ; done
cut -f1 x00_18.count.txt | sed '1d' | sort | uniq -d > duplicated.gene.list
sed '1iGeneid' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.list | csvtk join -j 4 -t -f '1;1' x00_18.count.txt - | cut -f1,7- > sheep.PCGlnc.gene.merged.count.txt
gene_num=`cut -f1 sheep.PCGlnc.gene.merged.count.txt | sed '1d' | wc -l`
sample_num=`head -1 sheep.PCGlnc.gene.merged.count.txt | awk '{print NF-1}'`
echo -e "#1.2\n${gene_num}\t${sample_num}" > sheep.PCGlnc.gene.merged.count.gct
sed 's/Geneid/Name/' sheep.PCGlnc.gene.merged.count.txt | awk '{ for (i=1; i<=NF; i++) { if (i==1) { printf $i"\t"$i"\t" } else if (i==NF) {printf $i"\n" } else { printf $i"\t" } } }' | sed 's/Name/Description/2' >> sheep.PCGlnc.gene.merged.count.gct
