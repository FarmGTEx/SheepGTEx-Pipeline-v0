import argparse
import pandas as pd

# --------------------------------------------------------
# Prepare QTL and matched SNP summary statistics
# - Load GWAS summary file
# - Extract SNPs for QTL and matched background
# - Format columns (A1/A2 uppercase, consistent IDs)
# - Output two files for LDAK input
# --------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Prepare QTL and matched SNP summary statistics for LDAK.")
    parser.add_argument("--gwas", required=True, help="GWAS summary statistics file")
    parser.add_argument("--qtl_snps", required=True, help="QTL SNP list file")
    parser.add_argument("--matched_snps", required=True, help="Matched SNP list file")
    parser.add_argument("--qtl_out", required=True, help="Output file for QTL SNPs")
    parser.add_argument("--matched_out", required=True, help="Output file for matched SNPs")
    args = parser.parse_args()

    # Load GWAS summary
    gwas = pd.read_csv(args.gwas, sep="\t")
    # Ensure column names are consistent
    gwas.columns = gwas.columns.str.strip()

    # Uppercase alleles
    gwas["A1"] = gwas["A1"].str.upper()
    gwas["A2"] = gwas["A2"].str.upper()

    # Load SNP lists
    qtl_snps = pd.read_csv(args.qtl_snps, sep="\t", header=None)[0].astype(str)
    matched_snps = pd.read_csv(args.matched_snps, sep="\t", header=None)[0].astype(str)

    # Extract subsets
    qtl_df = gwas[gwas["SNP"].isin(qtl_snps)].copy()
    matched_df = gwas[gwas["SNP"].isin(matched_snps)].copy()

    # Save outputs
    qtl_df.to_csv(args.qtl_out, sep="\t", index=False)
    matched_df.to_csv(args.matched_out, sep="\t", index=False)

if __name__ == "__main__":
    main()