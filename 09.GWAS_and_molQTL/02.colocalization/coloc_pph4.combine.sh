#!/bin/bash

# Input argument
tis=$1   # tissue name

# QTL directory
qtl_dir="QTL_results/${tis}"

# Loop over traits (skip log directory)
for trait in ${tis}/*; do
    if [ -d "$trait" ] && [ "$(basename $trait)" != "log" ]; then
        trait_name=$(basename $trait)
        echo "Processing tissue: $tis, trait: $trait_name"

        # Loop over QTL types
        for qtl in eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL; do
            pph4_file="${trait}/GWAS_${trait_name}_${qtl}.pph4"
            siglist_file="${tis}/${trait_name}.${qtl}.siglist"
            addcol_pph4_file="${trait}/GWAS_${trait_name}_${qtl}.addcol.pph4"
            addcol_sig_file="${trait}/GWAS_${trait_name}_${qtl}.addcol.sig.pph4"

            if [ -f "$pph4_file" ]; then
                echo "Found $pph4_file -> processing $qtl"

                # Add QTL and trait columns
                sed "s/^/${qtl}\t${trait_name}\t/g" "$pph4_file" | \
                    sed "1s/^$qtl\t$trait_name/qtl\ttrait/" > "$addcol_pph4_file"

                # Generate significant QTL list
                zcat ${qtl_dir}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.sig.txt.gz | \
                    cut -f1 | sed '1d' | sort -u > "$siglist_file"

                # Filter significant results
                if [ -s "$siglist_file" ]; then
                    sed '1iphenotype2' "$siglist_file" | \
                        csvtk join -t -f 'phenotype;phenotype2' "$addcol_pph4_file" - > "$addcol_sig_file"
                else
                    echo "Warning: $siglist_file is empty!"
                    rm -f "$addcol_sig_file"
                fi
            fi
        done
    fi
done

# Merge all traits × QTL results for this tissue
awk 'NR==1 || FNR>1' ${tis}/*/GWAS_*_*.addcol.pph4 > ${tis}/GWAS.addcol.pph4
awk 'NR==1 || FNR>1' ${tis}/*/GWAS_*_*.addcol.sig.pph4 > ${tis}/GWAS.addcol.sig.pph4
