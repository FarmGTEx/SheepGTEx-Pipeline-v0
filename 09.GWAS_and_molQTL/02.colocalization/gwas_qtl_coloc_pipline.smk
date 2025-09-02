# Config
TISSUE_LIST = "tissue40.list"       # tissue list file, first row = header
GWAS_BASE_DIR = "QTL_summaries"     # base dir for GWAS/QTL summaries
QTL_TYPES = ["eQTL", "sQTL", "eeQTL", "isoQTL", "stQTL", "3aQTL", "enQTL"]

# Load tissue names (skip header)
with open(TISSUE_LIST) as f:
    tissues = [line.strip().split()[0] for i, line in enumerate(f) if i > 0]

# Final outputs
rule all:
    input:
        expand("{tis}/GWAS.addcol.pph4", tis=tissues),
        expand("{tis}/GWAS.addcol.sig.pph4", tis=tissues),
        "GWAS.addcol.pph4",
        "GWAS.addcol.sig.pph4"

# Run coloc for each tissue × QTL
rule coloc_pph4:
    output:
        "{tis}/log/02.pph4.{tis}.{qtl}.log"
    params:
        gwas_dir=GWAS_BASE_DIR
    shell:
        """
        mkdir -p {wildcards.tis}/log
        bash coloc_pph4.sh {wildcards.tis} {wildcards.qtl} {params.gwas_dir} \
            > {output} 2>&1
        """

# Combine results per tissue
rule combine_tissue:
    input:
        expand("{tis}/log/02.pph4.{tis}.{qtl}.log", qtl=QTL_TYPES)
    output:
        "{tis}/GWAS.addcol.pph4",
        "{tis}/GWAS.addcol.sig.pph4"
    shell:
        """
        bash coloc_pph4.combine.sh {wildcards.tis}
        """

# Combine results across all tissues
rule combine_all:
    input:
        expand("{tis}/GWAS.addcol.pph4", tis=tissues),
        expand("{tis}/GWAS.addcol.sig.pph4", tis=tissues)
    output:
        "GWAS.addcol.pph4",
        "GWAS.addcol.sig.pph4"
    shell:
        """
        awk 'NR==1||FNR>1' */GWAS.addcol.pph4 > GWAS.addcol.pph4
        awk 'NR==1||FNR>1' */GWAS.addcol.sig.pph4 > GWAS.addcol.sig.pph4
        """
