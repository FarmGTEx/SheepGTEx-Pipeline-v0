### Usage: Expression quantification for all genes. Reference gtf file and sample bam files from STAR mapping are needed
rule gene_quant:
    input:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = GTFFILE
    output:
        smgtf = "04.MP1/gene/{sample}/{sample}.gtf",
        ae = "04.MP1/gene/{sample}/{sample}_stringtie.tsv",
        re = "04.MP1/gene/{sample}/{sample}_featureCounts.tsv"
    params:
        fqdir = "01.fastpQC"
    threads: 4
    shell:
        """
        stringtie -p {threads} -e -B -G {input.gtf} -o {output.smgtf} -A {output.ae} {input.bam}
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            featureCounts -T {threads} -p -t exon -g gene_id -a {input.gtf} -o {output.re} {input.bam}
        else
            featureCounts -T {threads} -t exon -g gene_id -a {input.gtf} -o {output.re} {input.bam}
        fi
        """

### Usage: Expression quantification for all exons, using featureCounts for read count quantification. Reference gtf file and sample bam files from STAR mapping are needed.
rule exon_quant:
    input:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = GTFFILE
    output:
        re = "04.MP1/exon/{sample}/{sample}_featureCounts.tsv"
    params:
        fqdir = "01.fastpQC"
    threads: 4
    shell:
        """
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            featureCounts -T {threads} -f -p -t exon -g gene_id -a {input.gtf} -o {output.re} {input.bam}
        else
            featureCounts -T {threads} -f -t exon -g gene_id -a {input.gtf} -o {output.re} {input.bam}
        fi
        """

#### Usage: Expression quantification for all enhancers, using featureCounts for read count quantification. Reference gtf file, enhancer saf file (please refer to featureCount Users Guide) and sample bam files from STAR mapping are needed.
rule enhancer_quant:
    input:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        enhancerSAF = ENHANCERSAF
    output:
        re = "04.MP1/enhancer/{sample}/{sample}_featureCounts.tsv"
    params:
        fqdir = "01.fastpQC"
    threads: 4
    shell:
        """
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            featureCounts -T {threads} -p -F SAF -a {input.enhancerSAF} -o {output.re} {input.bam}
        else
            featureCounts -T {threads} -F SAF -a {input.enhancerSAF} -o {output.re} {input.bam}
        fi
        """

### Usage: Expression quantification for transcripts/isoforms using salmon. Salmon can directly quantify transcript expression from fastq files,but remember to do the QC first. We need the gtf file and the genome fasta file to build the index for salmon first. Also we need the software salmon tools to generate the decoys file, which makes our result more accurate. More information about this checking this page https://github.com/COMBINE-lab/SalmonTools
rule transcript_quant:
    input:
        html = "01.fastpQC/{sample}/{sample}.html",
        json = "01.fastpQC/{sample}/{sample}.json",
        salmonidx = SALMONIDX
    output:
        sf = "04.MP1/transcript/{sample}/quant.sf",
        json = "04.MP1/transcript/{sample}/lib_format_counts.json"
    params:
        indir = "01.fastpQC",
        outdir = "04.MP1/transcript"
    threads: 4
    shell:
        """
        if [ -f "{params.indir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            salmon quant \
                -i {input.salmonidx} \
                -l A \
                -1 {params.indir}/{wildcards.sample}/{wildcards.sample}_1.clean.fq.gz \
                -2 {params.indir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz \
                -p {threads} \
                -o {params.outdir}/{wildcards.sample}
        else
            salmon quant \
                -i {input.salmonidx} \
                -l A \
                -r {params.indir}/{wildcards.sample}/{wildcards.sample}.clean.fq.gz \
                -p {threads} \
                -o {params.outdir}/{wildcards.sample}
        fi
        """

### Usage: Get the Bedgraph (.wig) format from bam, for 3' UTR polya estimation.
rule bam2Bedgraph:
    input:
        ibam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        ilog = "02.STAR2pass/{sample}/{sample}_Log.final.out"
    output:
        owig = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam.wig",
        ostat = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam.wig.mapped_reads.txt"
    threads: 1
    shell:
        """
        bedtools genomecov -bga -split -trackline -ibam {input.ibam} > {output.owig}
        depth=`grep 'Uniquely mapped reads number' {input.ilog} | cut -f2`
        echo -e "{wildcards.sample}\t{output.owig}\t$depth" > {output.ostat}
        """

rule combinewig:
    input:
        wigstats = expand("02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam.wig.mapped_reads.txt", sample=SAMPLES)
    output:
        allstats = "04.MP1/3UTRpolya/all.mapped_reads.txt"
    threads: 1
    shell:
        """
        cat {input.wigstats} > {output.allstats}
        """

### Usage: estimate 3'UTR polya using Dapars2 (https://github.com/3UTR/DaPars2/wiki). You need to convert the bam into Bedgraph file, and prepare config file first. (config file format: https://github.com/3UTR/DaPars2/wiki#step-3-run-dapars2-dapars2_multi_sample_multi_chrpy)
rule Dapars2_config:
    input:
        allstats = "04.MP1/3UTRpolya/all.mapped_reads.txt",
        threeUTRBED = THREEUTRBED
    output:
        stat = "04.MP1/3UTRpolya/all_mapped_reads.txt",
        config = "04.MP1/3UTRpolya/all.config"
    params:
        outdir = "04.MP1/3UTRpolya"
    threads: 10
    shell:
        """
        # generate sequencing depth file
        cut -f2,3 {input.allstats} > {output.stat}
        # generate config file
        wigs=`cut -f1 {output.stat} | tr '\\n' ',' | sed 's/,$/\\n/g'`
        echo -e "Annotated_3UTR={input.threeUTRBED}\\nAligned_Wig_files=$wigs\\nOutput_directory={params.outdir}/all\\nOutput_result_file=all\\nCoverage_threshold=10\\nNum_Threads={threads}\\nsequencing_depth_file={output.stat}" > {output.config}
        """

rule Dapars2:
    input:
        config = "04.MP1/3UTRpolya/all.config",
        dapars2DIR = DAPARS2DIR
    output:
        pdui = "04.MP1/3UTRpolya/all_{chrom}/all_result_temp.{chrom}.txt"
    threads: 10
    shell:
        """
        # run Dapars2
        python {input.dapars2DIR}/src/Dapars2_Multi_Sample.py {input.config} {wildcards.chrom}
        """

### Usage: estimate the RNA stability (https://github.com/daniel-munro/Pantry/blob/main/Project/steps/stability.smk). To do this, we should first extract exon and intron information from gtf file, and quantify them seperately. More information in https://github.com/csglab/CRIES.
rule consExon_quant:
    input:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        consExonGTF = CONSEXONGTF
    output:
        re = "04.MP1/RNA_stability/{sample}/{sample}.consExons.tsv"
    params:
        fqdir = "01.fastpQC"
    threads: 4
    shell:
        """
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            featureCounts -T {threads} -p -t exon -g gene_id -a {input.consExonGTF} --fracOverlap 1 -o {output.re} {input.bam}
        else
            featureCounts -T {threads} -t exon -g gene_id -a {input.consExonGTF} --fracOverlap 1 -o {output.re} {input.bam}
        fi
        """
        
rule intron_quant:
    input:
        bam = "02.STAR2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        intronGTF = INTRONGTF
    output:
        re = "04.MP1/RNA_stability/{sample}/{sample}.introns.tsv"
    params:
        fqdir = "01.fastpQC"
    threads: 4
    shell:
        """
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then 
            featureCounts -T {threads} -p -t intron -g gene_id -a {input.intronGTF} --fracOverlap 0 -o {output.re} {input.bam}
        else
            featureCounts -T {threads} -t intron -g gene_id -a {input.intronGTF} --fracOverlap 0 -o {output.re} {input.bam}
        fi
        """
