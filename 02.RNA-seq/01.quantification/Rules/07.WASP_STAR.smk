### Pipeline for mappability filtering using STAR (https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/src/run_STAR.py)
### (https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) (https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl/leafcutter#1-generating-wasp-corrected-alignments-with-star)
### (https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl/leafcutter#2-generating-exon-exon-junction-counts-with-regtools) and WASP (https://github.com/bmvdgeijn/WASP/tree/master/mapping)
rule wasp_star:
    input:
        gtf = GTFFILE,
        ibam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        vcf = "06.recalINFO/chrAuto.info.vcf.gz",
        staridx = STARIDX
    output:
        ovcf = temp("07.WASP_STAR/{sample}/{sample}.vcf"),
        obam = temp("07.WASP_STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"),
        logout = "07.WASP_STAR/{sample}/{sample}_Log.final.out",
        fbam = "07.WASP_STAR/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.bam"
    params:
        fqdir = "01.fastpQC",
        outdir = "07.WASP_STAR"
    threads: 4
    shell:
        """
        bcftools view {input.vcf} -s {wildcards.sample} --threads {threads} > {output.ovcf}
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then \
            # for pair-end reads
            STAR \
                --genomeDir {input.staridx} \
                --sjdbGTFfile {input.gtf} \
                --twopassMode Basic \
                --readFilesIn \
                    {params.fqdir}/{wildcards.sample}/{wildcards.sample}_1.clean.fq.gz \
                    {params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz \
                --outFileNamePrefix {params.outdir}/{wildcards.sample}/{wildcards.sample}_ \
                --runThreadN {threads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMismatchNmax 3 \
                --outSAMunmapped Within \
                --chimSegmentMin 10 \
                --chimOutType Junctions \
                --chimOutJunctionFormat 1 \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --varVCFfile {output.ovcf} \
                --waspOutputMode SAMtag \
                --outSAMattributes NH HI AS nM NM ch vA vG vW \
                --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample}
        else
            # for single-end reads
            STAR \
                --genomeDir {input.staridx} \
                --sjdbGTFfile {input.gtf} \
                --twopassMode Basic \
                --readFilesIn \
                    {params.fqdir}/{wildcards.sample}/{wildcards.sample}.clean.fq.gz \
                --outFileNamePrefix {params.outdir}/{wildcards.sample}/{wildcards.sample}_ \
                --runThreadN {threads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMismatchNmax 3 \
                --outSAMunmapped Within \
                --chimSegmentMin 10 \
                --chimOutType Junctions \
                --chimOutJunctionFormat 1 \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --varVCFfile {output.ovcf} \
                --waspOutputMode SAMtag \
                --outSAMattributes NH HI AS nM NM ch vA vG vW \
                --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample}
        fi
        # Filter out multi-mapping reads and reads that do not pass WASP filtering
        samtools view -h -q 255 {output.obam} | grep -v "vW:i:[2-7]" | samtools view -b > {output.fbam}
        """

### Generate junction file with Leafcutter (http://davidaknowles.github.io/leafcutter). Mapping bias should be removed first using WASP.
rule bam2junc:
    input:
        leafcutterDIR = LEAFCUTTERDIR,
        bam = "07.WASP_STAR/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.bam"
    output:
        junc = "07.WASP_STAR/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.bam.junc"
    threads: 4
    shell:
        """
        {input.leafcutterDIR}/scripts/bam2junc.sh {input.bam} {output.junc}
        """
