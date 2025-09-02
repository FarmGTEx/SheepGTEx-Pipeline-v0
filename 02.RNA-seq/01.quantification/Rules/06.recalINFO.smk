# recall the INFO field based on the final merged vcf. refnum: the number of samples in reference panel.
rule recalINFO:
    input:
        ibcf = "05.GLIMPSE2impute/{chrom}.bcf",
    output:
        ovcf = "06.recalINFO/{chrom}.vcf.gz"
    wildcard_constraints:
        chrom="chr\w+"
    params:
        scriptDIR = "Scripts"
    threads: 4
    shell:
        """
        python {params.scriptDIR}/recalINFO.py \
            --invcf {input.ibcf} | \
        bcftools +fill-tags -Oz -o {output.ovcf} --threads {threads} -- -t MAF
        bcftools index -f {output.ovcf} --threads {threads}
        """

rule info_filter:
    input:
        ivcf = "06.recalINFO/{chrom}.vcf.gz"
    output:
        ovcf = "06.recalINFO/{chrom}.info.vcf.gz"
    threads: 4
    shell:
        """
        bcftools view \
            -i "INFO>=0.75" \
            -O z -o {output.ovcf} {input.ivcf} --threads {threads}
        bcftools index -f {output.ovcf} --threads {threads}
        """

rule info_concat:
    input:
        ivcfs = expand("06.recalINFO/{chrom}.info.vcf.gz", chrom=CHROMS)
    output:
        autolist = "06.recalINFO/chrAuto.info.list",
        ovcf = "06.recalINFO/chrAuto.info.vcf.gz"
    threads: 4
    shell:
        """
        ls -1v {input.ivcfs} | grep -v 'chrX\|Y' > {output.autolist}
        bcftools concat -f {output.autolist} --threads {threads} -O z -o {output.ovcf}
        bcftools index {output.ovcf} --threads {threads}
        """

rule maf_filter:
    input:
        ivcf = "06.recalINFO/{chrom}.info.vcf.gz"
    output:
        ovcf = "06.recalINFO/{chrom}.filtered.vcf.gz"
    threads: 4
    shell:
        """
        bcftools view \
            -i "MAF>=0.05" \
            -O z -o {output.ovcf} {input.ivcf} --threads {threads}
        bcftools index -f {output.ovcf} --threads {threads}
        """

rule maf_concat:
    input:
        ivcfs = expand("06.recalINFO/{chrom}.filtered.vcf.gz", chrom=CHROMS)
    output:
        autolist = "06.recalINFO/chrAuto.filtered.list",
        ovcf = "06.recalINFO/chrAuto.filtered.vcf.gz"
    threads: 4
    shell:
        """
        ls -1v {input.ivcfs} | grep -v 'chrX\|Y' > {output.autolist}
        bcftools concat -f {output.autolist} --threads {threads} -O z -o {output.ovcf}
        bcftools index {output.ovcf} --threads {threads}
        """

