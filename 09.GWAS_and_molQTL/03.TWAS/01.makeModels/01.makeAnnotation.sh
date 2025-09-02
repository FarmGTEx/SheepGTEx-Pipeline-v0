#!/bin/bash
echo eQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p eQTL/${tis}/data eQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/covFile/${tis}.tsv eQTL/${tis}/data/covariates.txt ; done
csvtk join -t -H -f '4;1' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.bed /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.annot | sed 's/^chr//g' | awk -v OFS="\t" '{print $1,$4,$4,$2+1,$3,$5}' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > eQTL/gene_annot.txt

echo sQTL
# https://groups.google.com/g/predixcanmetaxcan/c/gin60-snZEM/m/Dyl4TSoBDQAJ
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p sQTL/${tis}/data sQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/covFile/${tis}.tsv sQTL/${tis}/data/covariates.txt ; done
csvtk join -t -H -f '1;2' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/01.MP1/02.stat/04.MP1/gene.list /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/03.MP2/splicing/All/All.leafcutter.clusters_to_genes.txt | sed 's/:/\t/g' | awk -v OFS="\t" '{print $2,$2":"$3":"$4":"$5":"$1,$2":"$3":"$4":"$5":"$1,$3,$4,"intron"}' | sed 's/^chr//g' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > sQTL/gene_annot.txt

echo eeQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p eeQTL/${tis}/data eeQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/covFile/${tis}.tsv eeQTL/${tis}/data/covariates.txt ; done
cut -f1 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/sheep.PCGlnc.exon.merged.count.txt | sed -e '1d' -e 's/:/\t/g' | awk -v OFS="\t" '{print $2,$1":"$2":"$3":"$4,$1":"$2":"$3":"$4,$3,$4,"exon"}' | sed 's/^chr//g' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > eeQTL/gene_annot.txt

echo isoQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p isoQTL/${tis}/data isoQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/covFile/${tis}.tsv isoQTL/${tis}/data/covariates.txt ; done
csvtk join -t -H -f '2;4' <(cut -f1 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/sheep.PCGlnc.isoform.merged.count.txt | sed -e '1d' -e 's/:/\t/g') <(awk '$3=="transcript"' /storage/public/home/2020060185/genome/sheep/reference/sheep.gtf | awk '{print $1"\t"$4"\t"$5"\t"$12}' | sed -e 's/"//g' -e 's/;//g') | awk -v OFS="\t" '{print $3,$1":"$2,$1":"$2,$4,$5,"isoform"}' | sed 's/^chr//g' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > isoQTL/gene_annot.txt

echo stQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p stQTL/${tis}/data stQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/covFile/${tis}.tsv stQTL/${tis}/data/covariates.txt ; done
cp eQTL/gene_annot.txt stQTL/gene_annot.txt

echo 3aQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p 3aQTL/${tis}/data 3aQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/covFile/${tis}.tsv 3aQTL/${tis}/data/covariates.txt ; done
cut -f1,4 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/03.MP2/3UTRpolya/Dapars2_res.all_chromosomes.txt | sed -e '1d' -e 's/|/\t/g' -e 's/:/\t/g' -e 's/\(.*\)-/\1\t/' | awk -v OFS="\t" '{print $3,$1"|"$2"|"$3"|"$4,$1"|"$2"|"$3"|"$4,$6,$7,"polya"}' | sed 's/^chr//g' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > 3aQTL/gene_annot.txt

echo enQTL
for tis in `cut -f1 tissue40.list | sed '1d'` ; do mkdir -p enQTL/${tis}/data enQTL/${tis}/output ; cp /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL_new/enQTL_new/${tis}/covFile/${tis}.tsv enQTL/${tis}/data/covariates.txt ; done
cut -f2 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL_new/enQTL_new/id.txt | sed -e '1d' -e 's/:/\t/g' -e 's/_/\t/' | awk  -v OFS="\t" '{print $3,$1":"$2"_"$3":"$4":"$5,$1":"$2"_"$3":"$4":"$5,$4,$5,"enhancer"}' | sed 's/^chr//g' | sed '1ichr\tgene_id\tgene_name\tstart\tend\tgene_type' > enQTL/gene_annot.txt
