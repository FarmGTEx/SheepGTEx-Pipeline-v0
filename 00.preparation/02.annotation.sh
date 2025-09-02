#!/bin/bash
## annovar annotation
gtfToGenePred sheep.gtf sheep_refGene.txt -genePredExt -allErrors
/storage/public/home/2021050427/software/annovar/retrieve_seq_from_fasta.pl sheep_refGene.txt --seqfile sheep.fna -format refGene -outfile sheep_refGeneMrna.fasta

## extract PCGlnc genes
awk '$3=="gene"' sheep.gtf > gene.gtf
grep 'gene_biotype "protein_coding";' gene.gtf  > pcg.gtf
grep 'gene_biotype "lncRNA";' gene.gtf  > lncRNA.gtf
cut -f9 pcg.gtf | sed 's/; /\t/g' | cut -f1 | sed 's/gene_id //g' | sed 's/"//g' > pcg.list
cut -f9 lncRNA.gtf | sed 's/; /\t/g' | cut -f1 | sed 's/gene_id //g' | sed 's/"//g' > lncRNA.list
cat pcg.gtf lncRNA.gtf > PCGlnc.gtf
awk '{print $1"\tprotein_coding"}' pcg.list > gene.annot
awk '{print $1"\tlncRNA"}' lncRNA.list >> gene.annot
cat pcg.list lncRNA.list > gene.list
awk '$3 =="gene" {gsub("\"","",$10);OFS="\t";print $1,int($4)-1,$5,$10,$7}' sheep.gtf | sed 's/;//g' | csvtk join -t -H -f '4;1' - gene.list > gene.bed

## extract PCGlnc genes with more than one exon
awk '$3=="exon"' sheep.gtf > exon.gtf
cut -f9 exon.gtf | sed 's/; /\t/g' | cut -f1 | sed -e 's/gene_id //g' -e 's/"//g' | csvtk join -t -H -f '1;1' gene.list - | uniq -c > gene_exon.count
awk '$1>1{print $2}' gene_exon.count > gene_exon2.list

## extract PCGlnc genes with more than one transcript
awk '$3=="transcript"' sheep.gtf > transcript.gtf
cut -f9 transcript.gtf | sed 's/; /\t/g' | sed -e 's/gene_id //g' -e 's/transcript_id //g' | sed 's/"//g' | cut -f1,2 | csvtk join -t -H -f '1;1' gene.list - > gene_transcript.list
sort -V gene_transcript.list | bedtools groupby -i - -g 1 -c 2 -o collapse | grep ',' | cut -f2 | sed 's/,/\n/g' > transcript2.list
csvtk join -t -H -f '2;1' gene_transcript.list transcript2.list > gene_transcript2.list