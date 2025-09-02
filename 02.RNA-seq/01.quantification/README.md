#### Pipeline and for the sheep&goat GTEx Project

### Software required (on linux):

## Miniconda3 (add mirrors to speed up: https://mirror.tuna.tsinghua.edu.cn/help/anaconda/)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

## Snakemake

conda install -n base -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal==7.30.1 pandas scipy click

## sratoolkit (prefetch)/ascp for data downloading

# For example, ascp:

wget https://download.asperasoft.com/download/sw/connect/3.10.0/ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz --no-check-certificate
tar zxvf ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.sh

## bedtools v2.31.0

wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
mv bedtools.static.binary bedtools ; chmod a+x bedtools ; mv bedtools ~/bin

## STAR v2.7.10b

wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b/bin/Linux_x86_64_static
ln -s $PWD/STAR ~/bin

## gatk v4.3.0.0

# Download gatk from the github releases page https://github.com/broadinstitute/gatk/releases

unzip gatk-4.3.0.0.zip
cd gatk-4.3.0.0
ln -s $PWD/gatk ~/bin

## fastp v0.23.2

wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
mv fastp ~/bin

## download GLIMPSE2 v2.0.0 (GLIMPSE2_*_static) from github release page https://github.com/odelaneau/glimpse/releases/

chmod a+x GLIMPSE2_*_static ; mv GLIMPSE2_*_static ~/bin

## stringtie v2.2.1

wget https://github.com/gpertea/stringtie/releases/download/v2.2.1/stringtie-2.2.1.Linux_x86_64.tar.gz
tar -xzf stringtie-2.2.1.Linux_x86_64.tar.gz
cd stringtie-2.2.1.Linux_x86_64
ln -s $PWD/stringtie ~/bin

## subread v2.0.5

wget https://udomain.dl.sourceforge.net/project/subread/subread-2.0.5/subread-2.0.5-Linux-x86_64.tar.gz --no-check-certificate
tar -xzf subread-2.0.5-Linux-x86_64.tar.gz
cd subread-2.0.5-Linux-x86_64/bin

# Add the software to the PATH in ~/.bashrc

export PATH=$PWD:$PATH

## bcftools, samtools and htslib v1.7+.

# According to http://www.htslib.org/download/

## Salmon v1.10.0

# Download salmon from github released page https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz

tar -xzf salmon-1.10.0_linux_x86_64.tar.gz
cd salmon-latest_linux_x86_64/bin
ln -s $PWD/salmon ~/bin

## SalmonTE v0.4 for TE quantification (https://github.com/hyunhwan-jeong/SalmonTE)
conda create -n SalmonTE pip
conda activate SalmonTE
pip3 install snakemake docopt pandas --user
conda deactivate

## gffread v0.12.7

# Download gffread from github released page https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz

tar -xzf gffread-0.12.7.Linux_x86_64.tar.gz
cd gffread-0.12.7.Linux_x86_64
ln -s $PWD/gffread ~/bin

## Download liftOver, gtfToGenePred, genePredToBed in https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64

wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod 755 liftOver gtfToGenePred genePredToBed
mv liftOver gtfToGenePred genePredToBed ~/bin

## DaPars2 for 3'UTR polya detection

# Downdload DaPars2 codes from github https://github.com/3UTR/DaPars2

unzip DaPars2-main.zip

## WASP for unbiased allele-specific read mapping

# NOTE: It's better to change the output mode into 'wb' in line 139-142 in mapping/find_intersecting_snps.py to output BAM file.

git clone https://github.com/bmvdgeijn/WASP.git
conda create -n WASP hdf5 pysam scipy
conda activate WASP
pip install tables -i https://pypi.tuna.tsinghua.edu.cn/simple
cd WASP/snp2h5 ; make
conda deactivate
ln -s $PWD/snp2h5 ~/bin
ln -s $PWD/fasta2h5 ~/bin

## LeafCutter v0.2.7 for splicing detection

# installation from http://davidaknowles.github.io/leafcutter/articles/Installation.html for splicing detection.

git clone https://github.com/davidaknowles/leafcutter

# install RegTools (https://regtools.readthedocs.io/en/latest/#installation)

git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
ln -s $PWD/regtools ~/bin

# install 

## phASER v1.1.1 for ASE (https://github.com/secastel/phaser)

git clone https://github.com/secastel/phaser.git
conda create -n phaser python=2.7 pysam scipy cython
conda activate phaser
cd phaser/phaser/
python setup.py build_ext --inplace
conda deactivate

### Indexed reference genome and gtf preparation, and also other file if needed (i.e. index of dbSNP). It will take several hours. Please modify the following if needed.

mkdir reference ; cd reference

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

## gatk v4.3.0.0 and samtools v1.7+ index

gatk CreateSequenceDictionary R=sheep.fna O=sheep.dict
samtools faidx sheep.fna
bedtools makewindows -g sheep.fna.fai -w 5000000 > sheep.5m.intervals
grep '^chr' sheep.5m.intervals | grep -v 'chrMT' > sheep.chr.5m.intervals

## salmon index

# get the cDNA sequence.

gffread sheep.gtf -g sheep.fna -w sheep.transcript.fa.tmp
cut -f 1 -d ' ' sheep.transcript.fa.tmp > sheep.transcript.fa ; rm sheep.transcript.fa.tmp

# The following command will output two files named gentrome.fa and decoys.txt, which will be used to build the index. The script can be downloaded from https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh. This will take ~50 minutes.

# To run the script, MashMap v2.0 is also needed (https://github.com/marbl/MashMap/releases/download/v2.0/mashmap-Linux64-v2.0.tar.gz).

bash ~/software/SalmonTools-master/scripts/generateDecoyTranscriptome.sh -a sheep.gtf -g sheep.fna -t sheep.transcript.fa -o salmonidx
salmon index -t salmonidx/gentrome.fa -d salmonidx/decoys.txt -i salmonidx/sheep.salmon

## Get the 3utr region for 3UTRpolya molecular phenotype analysis

gtfToGenePred -genePredExt sheep.gtf sheep.genepred
genePredToBed sheep.genepred sheep.genepred.bed
cut -f1,12 sheep.genepred > sheep.IDmapping.txt
python ~/software/DaPars2-master/src/DaPars_Extract_Anno.py -b sheep.genepred.bed -s sheep.IDmapping.txt -o sheep.3utr.bed

## Obtain the coordinates of intronic and constitutive exonic segments of the genes for mRNA stability detection (https://github.com/csglab/CRIES#step-1-creating-gtf-annotation-files). For the purpose of inferring stability from RNA-seq, we only consider constitutive exons, i.e. exons that are present in all isoforms of a unique gene. Codes edited from https://github.com/csglab/CRIES

grep -v '#' sheep.gtf | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[2]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\tncbi\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' > sheep.cons_exon.gtf
grep -v '#' sheep.gtf | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[2]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\tncbi\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' | bedtools sort -i stdin | awk -v FS='\t' '{ if( last_exon[$9]==1 && (last_exon_end[$9]+1)<($4-1) ) printf("%s\t%s\tintron\t%i\t%i\t%s\t%s\t%s\t%s\n",$1,$2,last_exon_end[$9]+1,$4-1,$6,$7,$8,$9); last_exon[$9]=1; last_exon_end[$9]=$5; }' > sheep.intron.gtf

## Prepare the enhancer saf file from cattle enhancer annotations.

enhancer_bed=bosTau9ToRamb2_E6.bed
ref_gtf=sheep.gtf
if [ ! -f enhancer.saf ]; then
grep -v '#' $ref_gtf | awk -v OFS="\t" '{print $1,$4,$5,$7}' | sort -k1,1 -k2,2n | uniq > ${ref_gtf}.sort.bed
    # Only include the enhancers in the gene body.
    bedtools intersect -a ${enhancer_bed} -b ${ref_gtf}.sort.bed -wa -wb | awk -v OFS="\t" '{print $1,$2,$3,$NF}' | uniq | awk '{{print "enhancer&"$1":"$2":"$3"\t"$1"\t"$2"\t"$3"\t"$4}}' | sed '1i GeneID\tChr\tStart\tEnd\tStrand' > enhancer.saf
fi

## Calculate the mappability to prepare the blacklist for ASE analysis. It will take several hours and it's better to submit the following steps into cluster by chromosome. For example in chr1:

chrom=chr1
# Step1. Separate the genome
samtools faidx ../sheep.fna ${chrom} > ${chrom}.fa
# Step2. Building the index
genmap index -F ${chrom}.fa -I ${chrom}_index
# Step3. Computing the mappability
genmap map -K 75 -E 2 -I ${chrom}_index -O ${chrom}_75_2 -t -w -bg
# Step4. Merge and extract the bed file with mappability<0.5
for chrom in chr{1..26} ; do cat ${chrom}_75_2.bedgraph ; done > chrAuto_75_2.bed
awk '{if ($4<0.5) print}' chrAuto_75_2.bed | bedtools sort -i - | bedtools merge -i - > chrAuto_75_2_0.5.bed

## Index of the dbSNP from reference panel

cd refpanel
dbSNP=Sheep3518.maf0.01.vcf
if [ ! -f ${dbSNP}.idx ]; then gatk IndexFeatureFile -I ${dbSNP}; fi

## WGS refpanel pre-process for GLIMPSE2 imputation. It will take several hours and it's better to submit the following steps into cluster by chromosome. For example in chr1:

# Step 0: Add the AN,AC tags if they are not exist. (~2 hours)

chr=chr1
bcftools +fill-tags split/Sheep3125.${chr}.BeaglePhase.rename.vcf.gz -Oz --threads 4 -o split/Sheep3125.${chr}.BeaglePhase.rename.AN_AC.vcf.gz -- -t AN,AC,MAF
bcftools index split/Sheep3125.${chr}.BeaglePhase.rename.AN_AC.vcf.gz

# Step 1: Chunk: create imputation chunks, for example in chr1 (~30 minutes).

refpanel=split/Sheep3125.${chr}.BeaglePhase.rename.AN_AC.vcf.gz
GLIMPSE2_chunk_static --input $refpanel --region ${chr} --window-mb 5.0 --buffer-mb 2.0 --output split/chunks.$chr.txt --sequential --threads 16

# Step 2: Split reference: create a binary reference panel for quick reading time. (~1 hours)

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    GLIMPSE2_split_reference_static --reference $refpanel --input-region ${IRG} --output-region ${ORG} --output ./split/split
done < split/chunks.$chr.txt

## Combine all the bins into a list file

for bin in `ls -1v split/split_chr*.bin` ; do basename $bin .bin ; done | sed 's/^split\_//g' > bins.txt
