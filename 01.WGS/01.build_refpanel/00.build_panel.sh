# 1. clean data (trimmomatic v0.39)

threads=$1
sampleid=$2
trimmomatic PE -threads ${threads} ${sampleid}_R1.fastq.gz ${sampleid}_R2.fastq.gz ${sampleid}_R1.paired.fastq.gz ${sampleid}_R1.unpaired.fastq.gz ${sampleid}_R2.paired.fastq.gz ${sampleid}_R2.unpaired.fastq.gz MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20

# 2. bwa align and sorted (bwa v0.7.5a-r405, samtools v1.9)

threads=$1
sampleid=$2
bwa mem -t ${threads} -R "@RG\tID:${sampleid}\tPL:illumina\tLB:${sampleid}\tSM:${sampleid}" ARS-UI_Ramb_v2.0_CM2CHR.fna ${sampleid}_R1.paired.fastq.gz ${sampleid}_R2.paired.fastq.gz | samtools view -S -b | samtools sort -@ ${threads} -m 4G -O bam -o ${sampleid}.bam

# 3. MarkDuplicates and generate index (GATK v4.1.4.1)

sampleid=$1
gatk MarkDuplicates -I ${sampleid}.bam	-M ${sampleid}.marked_dup_metrics.txt -O ${sampleid}_rmdup.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true

# 4. generate gvcf file for indi in chr region

sampleid=$1
chr=$2
gatk HaplotypeCaller --emit-ref-confidence GVCF -R ARS-UI_Ramb_v2.0_CM2CHR.fna -I ${sampleid}_rmdup.bam -L ${chr} -O ${sampleid}.${chr}.g.vcf.gz

# 5. CombineGVCFs for all samples and GenotypeGVCFs in single chr region
# Optional, if your sample size > 1000.

chr=$1
tmpdir=$2
threads=$3
gatk GenomicsDBImport -R ARS-UI_Ramb_v2.0_CM2CHR.fna --genomicsdb-workspace-path database_${chr} -V sample.${chr}.g.vcf.gz -L ${chr} --tmp-dir=${tmpdir} --reader-threads ${threads}

gatk GenotypeGVCFs -R ARS-UI_Ramb_v2.0_CM2CHR.fna -V gendb://database_${chr} -O Chr${chr}.vcf.gz

# 6. genotype in single chr and SelectVariants

chr=$1
gatk SelectVariants -R ARS-UI_Ramb_v2.0_CM2CHR.fna -V Chr${chr}.vcf.gz -select-type SNP -O Chr${chr}.raw.SNP.vcf.gz

# 7. mark the hard filter for snp

chr=$1
gatk VariantFiltration -R ARS-UI_Ramb_v2.0_CM2CHR.fna -V Chr${chr}.raw.SNP.vcf.gz -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName filtered -O Chr${chr}.SNP.mark.vcf.gz

# 8. select snps which passed the hard filter and was the binary snp for autosome (bcftools v1.9)

chr=$1
bcftools view Chr${chr}.SNP.mark.vcf.gz --min-af 0.01:minor -e 'F_MISSING>0.9' -m 2 -M 2 -v snps -f PASS -O z -o Chr${chr}.filter.maf0.01.missing.binary.vcf.gz

# 9. phase (beagle v5.4)

chr=$1
threads=$2
ne=$3
java -Xmx100g -jar beagle.22Jul22.46e.jar gt=Chr${chr}.filter.maf0.01.missing.binary.vcf.gz  nthreads={threads} ne=${ne} out=Chr${chr}.phased
