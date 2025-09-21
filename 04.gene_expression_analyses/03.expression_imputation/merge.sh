#!/bin/bash
# overlapped genes
for tis in `cut -f1 ../../tissue25.list | sed '1d'` ; do zcat ${tis}/phenotypes/${tis}.raw.expression.bed.gz | cut -f4 | sed '1d' ; done | sort | uniq -c > genecount.list
tissue_num=`sed '1d' ../../tissue25.list | wc -l`
awk -v tisc=$tissue_num '$1==tisc{print $2}' genecount.list | sed '1iphenotype_id' > overlapped.genelist
# merge all the tissues
for tis in `cut -f1 ../../tissue25.list | sed '1d' | head -1` ; do zcat ${tis}/phenotypes/${tis}.raw.expression.bed.gz | sed 's/^#//g' | csvtk join -j 4 -t -f 'phenotype_id;phenotype_id' - overlapped.genelist > express.tsv ; done
for tis in `cut -f1 ../../tissue25.list | sed '1,2d'` ; do zcat ${tis}/phenotypes/${tis}.raw.expression.bed.gz | sed 's/^#//g' | csvtk join -j 4 -t -f 'chr,start,end,phenotype_id;chr,start,end,phenotype_id' express.tsv - > tmp.tsv ; mv -f tmp.tsv express.tsv ; done
cut -f4- express.tsv | csvtk transpose -j 4 -d $'\t' - -D ',' | sed -e 's/phenotype_id/,tissue/' -e 's/_/,/' -e 's/_Embryo/,Embryo/g' -e 's/SAMN13567765,/SAMN13567765_/g' > express.csv
