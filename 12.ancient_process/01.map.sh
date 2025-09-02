AdapterRemoval --collapse --minadapteroverlap 1 --adapter1  \
--adapter2  --minlength 30 --gzip --trimns --trimqualities  \
	--file1 pair/${1}_R1.fastq.gz \
        --file2 pair/${1}_R2.fastq.gz    \
        --basename out/${1} 

cutadapt -O 1 -m 30 -a -o ${1}.s.fastq.gz Single-end/${1}.fastq.gz
bwa aln \
    reference.fasta \
    ${sample}.${library}.all.fastq.gz \
| bwa samse  \
    -r "@RG\tID:${sample}\tPL:illumina\tLB:${sample}\tSM:${sample}" \
    reference.fasta \
    - \
    ${sample}.${library}.all.fastq.gz > ${sample}.${library}.sam

samtools \
    sort \
    --reference reference.fasta \
    ${sample}.${library}.s.sam \
    --output-fmt BAM \
    -o ${sample}.${library}.bam
picard MarkDuplicates \
    I=${1}.merge.sort.bam \
    O=${1}.merge.sort.nodup.bam \
    M=${1}.marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT
samtools index ${1}.merge.sort.nodup.bam
java -jar GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R reference.fasta \
-I ${1}.merge.sort.nodup.bam \
-o relign_inter/${1}.intervals
java -jar GenomeAnalysisTK.jar \
   -T IndelRealigner \
       -R reference.fasta \
	-I ${1}.merge.sort.nodup.bam \
	-targetIntervals relign_inter/${1}.intervals \
	-o ${1}.relign_merge.bam 
mapDamage -i $BAM -r $REF -n 100000 -Q $BQ --no-stats --merge-reference-sequences -d ${FOLD}/${LIB}.mapDamage
samtools view -b -q 30 ${1}.relign_merge.bam -o ${1}.filtered.bam 
samtools index  ${1}.filtered.bam 


