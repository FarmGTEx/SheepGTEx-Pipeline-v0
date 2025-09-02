### Usage: Per-sample 2-pass mapping (STAR2passmode). STAR parameters are set according to STAR ENCODE standard options for long RNA-seq pipeline, human GTEx (https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md), cattle, pig and chicken GTEx.
### 
rule STAR2pass:
    input:
        html = "01.fastpQC/{sample}/{sample}.html",
        json = "01.fastpQC/{sample}/{sample}.json",
        gtf = GTFFILE,
        staridx = STARIDX
    output:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
        logout = "02.STAR2pass/{sample}/{sample}_Log.final.out",
        md5 = "02.STAR2pass/{sample}/md5.txt"
    params:
        indir = "01.fastpQC",
        outdir = "02.STAR2pass"
    threads: 8
    shell:
        """
        if [ -f "{params.indir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then \
            STAR \
                --genomeDir {input.staridx} \
                --sjdbGTFfile {input.gtf} \
                --twopassMode Basic \
                --readFilesIn \
                    {params.indir}/{wildcards.sample}/{wildcards.sample}_1.clean.fq.gz \
                    {params.indir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz \
                --outFileNamePrefix {params.outdir}/{wildcards.sample}/{wildcards.sample}_ \
                --runThreadN {threads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMismatchNmax 3 \ # this parameter should be removed or changed to 999 in the next version
                --outSAMunmapped Within \
                --chimSegmentMin 10 \
                --chimOutType Junctions \
                --chimOutJunctionFormat 1 \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outBAMcompression 10
        else
            STAR \
                --genomeDir {input.staridx} \
                --sjdbGTFfile {input.gtf} \
                --twopassMode Basic \
                --readFilesIn \
                    {params.indir}/{wildcards.sample}/{wildcards.sample}.clean.fq.gz \
                --outFileNamePrefix {params.outdir}/{wildcards.sample}/{wildcards.sample}_ \
                --runThreadN {threads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMismatchNmax 3 \ # this parameter should be removed or changed to 999 in the next version
                --outSAMunmapped Within \
                --chimSegmentMin 10 \
                --chimOutType Junctions \
                --chimOutJunctionFormat 1 \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outBAMcompression 10
        fi
        samtools index {output.bam} -@ {threads}
        md5sum {output.bam} > {output.md5}
        """
