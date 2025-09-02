### Usage: RNA-Seq imputation and phasing by population per bin using GLIMPSE2 (https://odelaneau.github.io/GLIMPSE/). The AN,AC tags should be in vcf files of refpanel (bcftools +fill-tags -Oz -o refpanl_AN.AC_vcf.gz -- refpanel.vcf.gz -t AN,AC). 
rule combinebam:
    input:
        bams = expand("03.mkdup/{sample}/{sample}_bqsr.cram", sample=SAMPLES)
    output:
        bamlist = "05.GLIMPSE2impute/all.bam.list"
    threads: 1
    shell:
        """
        for bam in {input.bams} ; do echo $bam `basename $bam _bqsr.cram` ; done > {output.bamlist}
        """

rule GLIMPSE2impute:
    input:
        bamlist = "05.GLIMPSE2impute/all.bam.list",
        ref = REFFA
    output:
        obcf = "05.GLIMPSE2impute/split/{dbbin}.bcf"
    params:
        bindir = BINDIR,
        binfile = "split_{dbbin}.bin"
    threads: 12
    shell:
        """
        GLIMPSE2_phase_static \
            --bam-list {input.bamlist} --reference {params.bindir}/{params.binfile} \
            --fasta {input.ref} --mapq 20 --ne 1000 \
            --output {output.obcf} --threads {threads}
        """

### Usage: Ligatation of multiple phased BCF/VCF files into a single whole chromosome file using GLIMPSE2 (https://odelaneau.github.io/GLIMPSE/). GLIMPSE2 is run in chunks that are ligated into chromosome-wide files maintaining the phasing.
rule GLIMPSE2ligate:
    input:
        ibcfs = expand("05.GLIMPSE2impute/split/{dbbin}.bcf", dbbin=DBBINS)
    output:
        filelist = "05.GLIMPSE2impute/{chrom}.txt",
        obcf = "05.GLIMPSE2impute/{chrom}.bcf"
    threads: 4
    shell:
        """
        ls -1v {input.ibcfs} | grep "{wildcards.chrom}_" > {output.filelist}
        GLIMPSE2_ligate_static \
            --input {output.filelist} \
            --output {output.obcf} \
            --threads {threads}
        """

