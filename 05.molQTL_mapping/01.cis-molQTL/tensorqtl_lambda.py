# -*- coding: utf-8 -*-
"""
# @FileName      : Untitled-1.py
# @Time          : 2024-03-14 00:30:41
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   :
"""

import click
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2
import tensorqtl
from tensorqtl import genotypeio, trans


def cal_lambdap(pvals):
    return np.nanmedian(norm.ppf(pvals/2)**2)/chi2.ppf(0.5, 1)

def cal_lambdaz(z):
    return np.nanmedian(z**2)/chi2.ppf(0.5, 1)


@click.command()
@click.option('--phenotype_bed_file')
@click.option('--covariates_file')
@click.option('--plink_prefix_path')
@click.option('--outfile')
@click.option('--gene_size', default=200, help='How many genes to input at once for QTL mapping')
@click.option('--batch_size', default=50000, help='chunk size for for GPU processing')
def main(phenotype_bed_file, covariates_file, plink_prefix_path, outfile, gene_size, batch_size):
    '''
    QTL mapping for all variant-phenotype pairs using tensorQTL and calculate the Inflation Factor (λ) for each gene
    '''
    # input phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates

    # input genotypes
    pr = genotypeio.PlinkReader(plink_prefix_path, select_samples=phenotype_df.columns)
    # load genotypes and variants into data frames
    genotype_df = pr.load_genotypes()

    # run and output
    with open(outfile, 'w') as fo:
        gene_num = phenotype_df.shape[0]
        for gene_start in range(0, gene_num, gene_size):
            gene_end = min(gene_start + gene_size, gene_num)
            print(f"Processing genes from {gene_start} to {gene_end}")

            # QTL-mapping (dense output)
            pval, b, b_se, af = trans.map_trans(genotype_df, phenotype_df.iloc[gene_start:gene_end, :], covariates_df,
                                                return_sparse=False, pval_threshold=1, maf_threshold=0, batch_size=batch_size)

            # output lambda
            for gene in phenotype_df.index[gene_start:gene_end]:
                lambdap = cal_lambdap(pval[gene])
                lambdaz = cal_lambdaz(b[gene]/b_se[gene])
                fo.write(f'{gene}\t{lambdap}\t{lambdaz}\n')

            del pval, b, b_se, af


if __name__ == '__main__':
    main()
