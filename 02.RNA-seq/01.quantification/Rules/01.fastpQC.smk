### Usage: Quality control of merged raw reads (per sample) using fastp.
### 
rule fastpQC:
    input:
        md5 = "00.mergefq/{sample}/md5.txt"
    output:
        html = "01.fastpQC/{sample}/{sample}.html",
        json = "01.fastpQC/{sample}/{sample}.json",
        md5o = "01.fastpQC/{sample}/md5.txt"
    params:
        indir = "00.mergefq",
        outdir = "01.fastpQC"
    threads: 4
    shell:
        """
        if [ -f "{params.indir}/{wildcards.sample}/{wildcards.sample}_2.fq.gz" ];then \
            fastp \
                -i {params.indir}/{wildcards.sample}/{wildcards.sample}_1.fq.gz \
                -o {params.outdir}/{wildcards.sample}/{wildcards.sample}_1.clean.fq.gz \
                -I {params.indir}/{wildcards.sample}/{wildcards.sample}_2.fq.gz \
                -O {params.outdir}/{wildcards.sample}/{wildcards.sample}_2.clean.fq.gz \
                -f 3 -t 3 -l 36 -r -W 4 -M 15 \
                --detect_adapter_for_pe \
                --html {output.html} \
                --json {output.json} \
                --thread {threads} --compression 9
        else 
            fastp \
                -i {params.indir}/{wildcards.sample}/{wildcards.sample}.fq.gz \
                -o {params.outdir}/{wildcards.sample}/{wildcards.sample}.clean.fq.gz \
                -f 3 -t 3 -l 36 -r -W 4 -M 15 \
                --html {output.html} \
                --json {output.json} \
                --thread {threads} --compression 9
        fi
        md5sum {params.outdir}/{wildcards.sample}/{wildcards.sample}*.clean.fq.gz > {output.md5o}
        """

