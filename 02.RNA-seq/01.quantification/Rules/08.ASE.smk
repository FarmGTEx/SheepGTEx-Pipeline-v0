### Generate ASE Data with phASER (https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser). Mapping bias should be removed first using WASP.
rule phASER:
    input:
        wasppy = WASPPY,
        waspDIR = WASPDIR,
        phaserDIR = PHASERDIR,
        blacklist = BLACKLIST,
        fbam = "07.WASP_STAR/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.bam",
        vcf = "06.recalINFO/chrAuto.filtered.vcf.gz"
    output:
        fdbam = temp("08.ASE/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.rmdup.bam"),
        fdsbam = temp("08.ASE/{sample}/{sample}_Aligned.sortedByCoord.out.filtered.rmdup.sorted.bam"),
        ovcf = temp("08.ASE/{sample}/{sample}.filtered.vcf.gz"),
        ha_counts = "08.ASE/{sample}/{sample}.haplotypic_counts.txt"
    threads: 4
    params:
        fqdir = "01.fastpQC",
        outdir = "08.ASE"
    shell:
        """
        bcftools view {input.vcf} -s {wildcards.sample} -Oz --threads {threads} -o {output.ovcf}
        bcftools index -t {output.ovcf} --threads {threads}
        if [ -f "{params.fqdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz" ];then \
            # Filter duplicate reads for pair-end reads
            {input.wasppy} {input.waspDIR}/mapping/rmdup_pe.py \
                {input.fbam} {output.fdbam}
            samtools sort -o {output.fdsbam} --threads {threads} {output.fdbam}
            samtools index {output.fdsbam} -@ {threads}
            # phASER
            source /storage/public/home/2020060185/anaconda3/envs/phaser/bin/activate phaser
            python {input.phaserDIR}/phaser/phaser.py \
                --bam {output.fdsbam} \
                --vcf {output.ovcf} \
                --paired_end 1 \
                --mapq 255 --baseq 10 \
                --sample {wildcards.sample} \
                --haplo_count_blacklist {input.blacklist} \
                --threads {threads} \
                --o {params.outdir}/{wildcards.sample}/{wildcards.sample} \
                --output_read_ids 0 \
                --gw_phase_vcf 1 \
                --pass_only 0
        else
            # Filter duplicate reads for single-end reads
            {input.wasppy} {input.waspDIR}/mapping/rmdup.py \
                {input.fbam} {output.fdbam}
            samtools sort -o {output.fdsbam} --threads {threads} {output.fdbam}
            samtools index {output.fdsbam} -@ {threads}
            # phASER
            source /storage/public/home/2020060185/anaconda3/envs/phaser/bin/activate phaser
            python {input.phaserDIR}/phaser/phaser.py \
                --bam {output.fdsbam} \
                --vcf {output.ovcf} \
                --paired_end 0 \
                --mapq 255 --baseq 10 \
                --sample {wildcards.sample} \
                --haplo_count_blacklist {input.blacklist} \
                --threads {threads} \
                --o {params.outdir}/{wildcards.sample}/{wildcards.sample} \
                --output_read_ids 0 \
                --gw_phase_vcf 1 \
                --pass_only 0
        fi
        """

# Produce gene level haplotype counts for allelic expression studies
rule phASER_gene_ae:
    input:
        phaserDIR = PHASERDIR,
        genebed = GENEBED,
        ha_counts = "08.ASE/{sample}/{sample}.haplotypic_counts.txt"
    output:
        gene_ae = "08.ASE/{sample}/{sample}_phaser.gene_ae.txt"
    threads: 1
    shell:
        """
        source /storage/public/home/2020060185/anaconda3/envs/phaser/bin/activate phaser
        python {input.phaserDIR}/phaser_gene_ae/phaser_gene_ae.py \
            --haplotypic_counts {input.ha_counts} --features {input.genebed} \
            --o {output.gene_ae}
        """
