# -*- coding: utf-8 -*-
"""
# @FileName      : mediation.py
# @Time          : 2026-3-19 15:09:51
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   : trans-eQTL mediation analysis of cis-eQTL within tissue.
"""

import click
import pandas as pd
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
import pysam

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# extract cis-SNPs
def load_cis_genotype(vcf, chrom, start, end, samples):    
    geno_list = []
    snp_ids = []

    for record in vcf.fetch(chrom, start, end):
        genotypes = []
        for s in samples:
            gt = record.samples[s]["GT"]
            # 转成 0/1/2
            if None in gt:
                g = np.nan
            else:
                g = sum(gt)
            genotypes.append(g)
        
        geno_list.append(genotypes)
        snp_ids.append(record.id)

    Z = pd.DataFrame(geno_list, columns=samples, index=snp_ids).T
    return Z  # samples × SNPs

# orthogonalization
def residualize(y, covariates):
    X = sm.add_constant(covariates)
    model = sm.OLS(y, X).fit()
    return model.resid


# tsls model
def tsls_mediation(Z, x, y, alpha=100):
    """
    Z: (n_samples, n_snps) cis窗口基因型矩阵
    x: (n_samples,) cis gene expression
    y: (n_samples,) trans gene expression
    """

    # ---------- Step 1: Ridge regression ----------
    ridge = Ridge(alpha=alpha)
    ridge.fit(Z, x)

    # predicted cis expression
    x_hat = ridge.predict(Z)

    # ---------- Step 2: TSLS regression ----------
    X = sm.add_constant(x_hat)
    model = sm.OLS(y, X).fit()

    beta_tsls = model.params[1]
    pval = model.pvalues[1]

    return beta_tsls, pval, x_hat


# permutation
def tsls_with_permutation(Z, x, y, n_perm=100, alpha=100):
    beta_real, p_real, x_hat = tsls_mediation(Z, x, y, alpha)

    perm_pvals = []
    for _ in range(n_perm):
        y_perm = np.random.permutation(y)
        X = sm.add_constant(x_hat)
        model = sm.OLS(y_perm, X).fit()
        perm_pvals.append(model.pvalues[1])

    # empirical p-value
    emp_p = np.mean(np.array(perm_pvals) <= p_real)

    return beta_real, p_real, emp_p


@click.command()
@click.option('--info', help="information file of lead trans-eQTL and cis-eGene, including following fields: tissue, trans_gene, variant_id, cis_gene")
@click.option('--tissue', help="Tissue name")
@click.option('--cis_expr', help="Gene expression file used for cis-eQTL mapping")
@click.option('--trans_expr', help="Gene expression file used for trans-eQTL mapping")
@click.option('--vcf', help="Genotype file used for cis-eQTL mapping")
@click.option('--covar', help="Covariate file used for eQTL mapping")
@click.option('--out', help="out file")
def main(info, tissue, cis_expr, trans_expr, vcf, covar, out):
    # input data
    dfi = pd.read_csv(info, sep='\t')
    tis = tissue
    cis_expr = pd.read_csv(cis_expr, sep='\t', index_col=3)
    trans_expr = pd.read_csv(trans_expr, sep='\t', index_col=3)
    vcf = pysam.VariantFile(vcf)
    covar = pd.read_csv(covar, sep="\t", index_col=0)
    samples = trans_expr.columns[3:]

    # mediation analysis for each lead trans-eQTL within tissue
    results = []
    for i, row in dfi[dfi['tissue']==tis].iterrows():
        variant_id = row.variant_id
        trans_gene = row.trans_gene
        cis_gene = row.cis_gene
        chrom = cis_expr.loc[cis_gene, '#chr']
        start = cis_expr.loc[cis_gene, 'end'] - 1000000
        end = cis_expr.loc[cis_gene, 'end'] + 1000000
        
        # cis-SNPs ±1Mb
        Z = load_cis_genotype(vcf, chrom, start, end, samples)

        # expression
        x = pd.to_numeric(cis_expr.loc[cis_gene, samples])
        y = pd.to_numeric(trans_expr.loc[trans_gene, samples])

        # covariates
        cov_sub = covar.T.loc[samples]

        # orthogonalization
        x_res = residualize(x, cov_sub)
        y_res = residualize(y, cov_sub)

        # TSLS mediation
        beta, pval, emp_p = tsls_with_permutation(Z, x_res, y_res)
        results.append([tis, trans_gene, variant_id, cis_gene, beta, pval, emp_p])

    # save data
    results_df = pd.DataFrame(results, columns=['tissue', 'trans_gene', 'variant_id', 'cis_gene', 'beta', 'pval', 'emp_p'])
    rej, qvals  = fdrcorrection(results_df['emp_p'].values, alpha=0.05)
    results_df['qvalue'] = qvals
    results_df['significant'] = rej
    results_df.to_csv(out, index=False)

if __name__ == '__main__':
    main()
