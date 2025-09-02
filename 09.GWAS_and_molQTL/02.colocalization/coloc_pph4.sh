#!/bin/bash

# Arguments
tis=$1        # tissue name
qtl=$2        # QTL type
gwas_base_dir=$3  # base directory of GWAS results

# QTL directory
qtl_dir="QTL_results/${tis}"

# Loop over GWAS traits
for trait_dir in ${gwas_base_dir}/*; do
    [ -d "$trait_dir" ] || continue
    trait=$(basename "$trait_dir")
    gwas_dir="${trait_dir}"

    echo "[Processing] $tis - $qtl - $trait"

    # Output directory and file
    out_dir="${tis}/${trait}"
    mkdir -p "$out_dir"
    output_file="${out_dir}/GWAS_${trait}_${qtl}.pph4"
    echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tSNP\tphenotype\tchrom\ttissue" > "$output_file"

    # Loop over chromosomes
    for chr in {1..26}; do
        gwas_file="${gwas_dir}/chr${chr}.tsv"
        qtl_file="${qtl_dir}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.chr${chr}.txt.gz"
        output_prefix="${out_dir}/GWAS_${trait}_${qtl}.pph4"

        # Skip if GWAS file does not exist
        [ -f "$gwas_file" ] || continue

        # Run R script
        Rscript scripts/gwas_coloc_pph4.R \
            "$gwas_file" "$qtl_file" "$tis" "chr${chr}" "$output_prefix"
    done
done
