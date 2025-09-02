### Usage: Bam AddOrReplaceReadGroups, MarkDuplicates, SplitNCigarReads and BaseRecalibrator, using gatk (v4.3.0.0).
rule AddOrReplaceReadGroups:
    input:
        ibam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        obam = temp("03.mkdup/{sample}/{sample}_rg.bam")
    threads: 4
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads}" AddOrReplaceReadGroups \
            --spark-runner LOCAL \
            -I {input.ibam} -O {output.obam} \
            -LB {wildcards.sample} -PL ILLUMINA -PU {wildcards.sample} \
            -SM {wildcards.sample} -ID {wildcards.sample} -FO {wildcards.sample} \
            --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
        """

rule MarkDuplicates:
    input:
        ibam = "03.mkdup/{sample}/{sample}_rg.bam"
    output:
        obam = temp("03.mkdup/{sample}/{sample}_mkdup.bam"),
        metrics = "03.mkdup/{sample}/{sample}_metrics.txt"
    threads: 4
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads} -XX:-UseGCOverheadLimit" MarkDuplicates \
            --spark-runner LOCAL \
            -I {input.ibam} -O {output.obam} -M {output.metrics} \
            --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
            --REMOVE_DUPLICATES true
        """

rule SplitNCigarReads:
    input:
        ibam = "03.mkdup/{sample}/{sample}_mkdup.bam",
        ref = REFFA
    output:
        obam = temp("03.mkdup/{sample}/{sample}_cigar.bam")
    threads: 4
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads}" SplitNCigarReads \
            --spark-runner LOCAL \
            -I {input.ibam} -O {output.obam} -R {input.ref} \
            --create-output-bam-index true
        """

rule BaseRecalibrator:
    input:
        ibam = "03.mkdup/{sample}/{sample}_cigar.bam",
        ref = REFFA,
        dbsnp = DBSNP  # vcf file should be indexed by gatk
    output:
        bqsr = "03.mkdup/{sample}/{sample}_bqsr.table"
    threads: 4
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads}" BaseRecalibrator \
            --spark-runner LOCAL \
            -I {input.ibam} -R {input.ref} -O {output.bqsr} \
            --known-sites {input.dbsnp}
        """

rule ApplyBQSR:
    input:
        ibam = "03.mkdup/{sample}/{sample}_cigar.bam",
        bqsr = "03.mkdup/{sample}/{sample}_bqsr.table"
    output:
        obam = temp("03.mkdup/{sample}/{sample}_bqsr.bam")
    threads: 4
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads}" ApplyBQSR \
            --spark-runner LOCAL \
            -I {input.ibam} --bqsr-recal-file {input.bqsr} \
            -O {output.obam} \
            --create-output-bam-index true --create-output-bam-md5 true
        """

rule bqsrbam2cram:
    input:
        ibam = "03.mkdup/{sample}/{sample}_bqsr.bam",
        ref = REFFA
    output:
        ocram = "03.mkdup/{sample}/{sample}_bqsr.cram",
        md5 = "03.mkdup/{sample}/md5.txt"
    threads: 4
    shell:
        """
        samtools view {input.ibam} -T {input.ref} -C -o {output.ocram} \
            --write-index --threads {threads}
        md5sum {output.ocram} > {output.md5}
        """
