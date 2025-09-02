#!/bin/bash
## Download the reference data using ascp. Here is an example for the preparation of sheep reference genome and gtf file

ascp -v -k 1 -T -l 500m -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/genomes/all/GCF/016/772/045/GCF_016772045.1_ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna.gz .
ascp -v -k 1 -T -l 500m -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/genomes/all/GCF/016/772/045/GCF_016772045.1_ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz .

# Unzip

gunzip -c GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna.gz > GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna
gunzip -c GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz > GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf

## Change the chromosome name

echo '''{
s/NC_056054.1/chr1/
s/NC_056055.1/chr2/
s/NC_056056.1/chr3/
s/NC_056057.1/chr4/
s/NC_056058.1/chr5/
s/NC_056059.1/chr6/
s/NC_056060.1/chr7/
s/NC_056061.1/chr8/
s/NC_056062.1/chr9/
s/NC_056063.1/chr10/
s/NC_056064.1/chr11/
s/NC_056065.1/chr12/
s/NC_056066.1/chr13/
s/NC_056067.1/chr14/
s/NC_056068.1/chr15/
s/NC_056069.1/chr16/
s/NC_056070.1/chr17/
s/NC_056071.1/chr18/
s/NC_056072.1/chr19/
s/NC_056073.1/chr20/
s/NC_056074.1/chr21/
s/NC_056075.1/chr22/
s/NC_056076.1/chr23/
s/NC_056077.1/chr24/
s/NC_056078.1/chr25/
s/NC_056079.1/chr26/
s/NC_056080.1/chrX/
s/NC_001941.1/chrMT/
}''' > rename.sed
sed -f rename.sed GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna > sheep.fna # It will take ~5 minutes.
sed -f rename.sed GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf | sed 's/transcript_id ""; //g' > sheep.gtf

## add Y chromosome (CM022046.1) to the reference

cat husheep_Y.fa >> sheep.fna
python chrY_gff2gtf.py --ingff husheep_Y.gff --outgtf husheep_Y.gtf
cat husheep_Y.gtf >> sheep.gtf

## build reference index using STAR. It needs at least 40G memory and will take ~20 minutes.

STAR --runMode genomeGenerate --genomeDir staridx --runThreadN 4 --genomeFastaFiles sheep.fna --sjdbGTFfile sheep.gtf

## BWA index

bwa index sheep.fna

## BOWTIE (v1.3.1) index

bowtie-build -f sheep.fna --threads 16 sheep

## gatk v4.3.0.0 and samtools v1.7+ index

gatk CreateSequenceDictionary -R sheep.fna -O sheep.dict
samtools faidx sheep.fna
bedtools makewindows -g sheep.fna.fai -w 5000000 | awk '{print $1":"$2+1"-"$3}' > sheep.5m.intervals
grep '^chr' sheep.5m.intervals | grep -v 'chrMT' > sheep.chr.5m.intervals

## salmon index

# get the cDNA sequence.

gffread sheep.gtf -g sheep.fna -w sheep.transcript.fa.tmp
cut -f 1 -d ' ' sheep.transcript.fa.tmp > sheep.transcript.fa ; rm sheep.transcript.fa.tmp

# The following command will output two files named gentrome.fa and decoys.txt, which will be used to build the index. The script can be downloaded from https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh. This will take ~50 minutes.

# To run the script, MashMap v2.0 is also needed (https://github.com/marbl/MashMap/releases/download/v2.0/mashmap-Linux64-v2.0.tar.gz).

bash ~/software/SalmonTools-master/scripts/generateDecoyTranscriptome.sh -J 32 -a sheep.gtf -g sheep.fna -t sheep.transcript.fa -o salmonidx
salmon index -p 32 -t salmonidx/gentrome.fa -d salmonidx/decoys.txt -i salmonidx/sheep.salmon

## Get the 3utr region for 3UTRpolya molecular phenotype analysis

gtfToGenePred -genePredExt sheep.gtf sheep.genepred
genePredToBed sheep.genepred sheep.genepred.bed
cut -f1,12 sheep.genepred > sheep.IDmapping.txt
python ~/software/DaPars2-master/src/DaPars_Extract_Anno.py -b sheep.genepred.bed -s sheep.IDmapping.txt -o sheep.3utr.bed

## Obtain the coordinates of intronic and constitutive exonic segments of the genes for mRNA stability detection (https://github.com/csglab/CRIES#step-1-creating-gtf-annotation-files). For the purpose of inferring stability from RNA-seq, we only consider constitutive exons, i.e. exons that are present in all isoforms of a unique gene. Codes edited from https://github.com/csglab/CRIES

grep -v '#' sheep.gtf | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[2]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\tncbi\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' > sheep.cons_exon.gtf
grep -v '#' sheep.gtf | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[2]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\tncbi\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' | bedtools sort -i stdin | awk -v FS='\t' '{ if( last_exon[$9]==1 && (last_exon_end[$9]+1)<($4-1) ) printf("%s\t%s\tintron\t%i\t%i\t%s\t%s\t%s\t%s\n",$1,$2,last_exon_end[$9]+1,$4-1,$6,$7,$8,$9); last_exon[$9]=1; last_exon_end[$9]=$5; }' > sheep.intron.gtf

## Prepare the enhancer saf file from sheep enhancer annotations.

enhancer_bed=E6_Gs.bed
ref_gtf=sheep.gtf
if [ ! -f enhancer.saf ]; then
grep -v '#' $ref_gtf | awk -v OFS="\t" '{print $1,$4,$5,$7}' | sort -k1,1 -k2,2n | uniq > ${ref_gtf}.sort.bed
    # Only include the enhancers outside the gene body.
    bedtools intersect -a ${enhancer_bed} -b ${ref_gtf}.sort.bed -wa -v | awk -v OFS="\t" '{print $1,$2,$3}' | uniq | awk '{{print "enhancer_"$1":"$2":"$3"\t"$1"\t"$2"\t"$3"\t."}}' | sed '1i GeneID\tChr\tStart\tEnd\tStrand' > enhancer.saf
fi

