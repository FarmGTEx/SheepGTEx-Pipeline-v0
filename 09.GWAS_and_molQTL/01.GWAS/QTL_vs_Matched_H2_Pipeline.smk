# ===============================================
# Snakemake pipeline for LDAK heritability analysis
# Step1: SNP tagging per chromosome
# Step2: Prepare GWAS summary stats (QTL / matched)
# Step3: Heritability estimation (QTL / matched)
# ===============================================

CHRS = [str(i) for i in range(1, 27)]  # chromosomes 1-26
TRAITS = ["trait1", "trait2"]  # example, replace with actual traits

# ================
# Final target
# ================
rule all:
    input:
        expand("results/heritability/{trait}/qtl_h2.txt", trait=TRAITS),
        expand("results/heritability/{trait}/matched_h2.txt", trait=TRAITS)

# ================================
# Step1: Tagging calculation
# ================================
rule tagging:
    input:
        bfile = "data/genotype/chr{chr}"  # plink bed/bim/fam prefix
    output:
        tag = "results/tagging/chr{chr}_tagging.tagging"
    threads: 4
    shell:
        """
        echo "[INFO] Calculating tagging for chr{wildcards.chr}"
        ldak5.2 --calc-tagging {output.tag} \
                --bfile {input.bfile} \
                --power -0.25 \
                --window-kb 1000
        """

# ================================
# Step2: Prepare GWAS summary statistics
# ================================
rule prepare_sumstats:
    input:
        gwas = "data/gwas/{trait}.txt",
        tagging = expand("results/tagging/chr{chr}_tagging.tagging", chr=CHRS)
    output:
        qtl_sum = "results/sumstats/{trait}_qtl.sumstats",
        matched_sum = "results/sumstats/{trait}_matched.sumstats"
    shell:
        """
        mkdir -p results/sumstats/{wildcards.trait}

        # Extract QTL SNPs and matched SNPs
        python scripts/prepare_qtl_matched_sumstats.py \
            --gwas {input.gwas} \
            --tagging results/tagging/ \
            --qtl_out {output.qtl_sum} \
            --matched_out {output.matched_sum}
        """

# ================================
# Step3a: Heritability for QTL SNPs
# ================================
rule heritability_qtl:
    input:
        sumstats = "results/sumstats/{trait}_qtl.sumstats",
        tagging = expand("results/tagging/chr{chr}_tagging.tagging", chr=CHRS)
    output:
        "results/heritability/{trait}/qtl_h2.txt"
    shell:
        """
        mkdir -p results/heritability/{wildcards.trait}
        echo "[INFO] Running LDAK for QTL SNPs ({wildcards.trait})"
        ldak5.2 --sum-hers {output} \
                --tagfile results/tagging/chr*_tagging.tagging \
                --summary {input.sumstats} \
                --cutoff 0.01 \
                --check-sums NO
        """

# ================================
# Step3b: Heritability for matched SNPs
# ================================
rule heritability_matched:
    input:
        sumstats = "results/sumstats/{trait}_matched.sumstats",
        tagging = expand("results/tagging/chr{chr}_tagging.tagging", chr=CHRS)
    output:
        "results/heritability/{trait}/matched_h2.txt"
    shell:
        """
        mkdir -p results/heritability/{wildcards.trait}
        echo "[INFO] Running LDAK for matched SNPs ({wildcards.trait})"
        ldak5.2 --sum-hers {output} \
                --tagfile results/tagging/chr*_tagging.tagging \
                --summary {input.sumstats} \
                --cutoff 0.01 \
                --check-sums NO
        """
