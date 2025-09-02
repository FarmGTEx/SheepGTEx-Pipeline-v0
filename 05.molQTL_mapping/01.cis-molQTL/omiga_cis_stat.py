# -*- coding: utf-8 -*-
"""
# @FileName      : omiga_cis_stat.py
# @Time          : 2024-06-20 18:08:29
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   :
"""

import click
import numpy as np
import pandas as pd
import scipy.stats as stats
import glob


@click.command()
@click.option('--chrlist', help='list of chromosomes to be analysed')
@click.option('--indir', help='input directory of tensorqtl and omiga results')
@click.option('--outfile', help='output file')
def main(chrlist, indir, outfile):
    '''
    Calculate the correlation between cis-QTLs of LM (tensorqtl) and LMM (OmiGA) for each gene.
    '''
    tissue = indir.strip().split('/')[0]
    permfiles = glob.glob(f'{indir}/tensorqtl/permutation/*.cis_qtl_fdr0.05.txt.gz')
    if len(permfiles)!=1:
        raise IOError(f"Need one permutation file, but {permfiles} detected.")
    perm = pd.read_csv(permfiles[0], sep='\t')
    perm['z'] = perm['slope']/perm['slope_se']

    with open(outfile, 'w') as f:
        for chrom in [line.strip() for line in open(chrlist, 'r')]:
            lmfiles = glob.glob(f'{indir}/tensorqtl/nominal/*.cis_qtl_pairs.{chrom}.txt.gz')
            if len(lmfiles)!=1:
                raise IOError(f"Need one nominal file, but {lmfiles} detected.")
            dflm = pd.read_csv(lmfiles[0], sep='\t')
            dflm['z'] = dflm['slope']/dflm['slope_se']
            genelist = dflm['phenotype_id'].unique()

            lmmfiles = glob.glob(f'{indir}/omiga/cis_lmm/*.cis_qtl_pairs.{chrom}.txt.gz')
            if len(lmmfiles)!=1:
                raise IOError(f"Need one nominal file, but {lmmfiles} detected.")
            dflmm = pd.read_csv(lmmfiles[0], sep='\t')
            dflmm['Z'] = dflmm['beta_g1']/dflmm['beta_se_g1']

            for gene in genelist:

                ## calculate correlation of nominal results between LM and LMM
                df = pd.merge(dflm[dflm['phenotype_id']==gene].dropna(), dflmm[dflmm['pheno_id']==gene].dropna(), left_on=['phenotype_id','variant_id'], right_on=['pheno_id', 'variant_id'])
                rs = rz = rp = np.nan
                cis_num = df.shape[0]
                if df.shape[0]>1:
                    rs, Ps = stats.pearsonr(df['slope'], df['beta_g1'])
                    rz, Pz = stats.pearsonr(df['z'], df['Z'])
                    rp, Pp = stats.pearsonr(-np.log10(df['pval_nominal']), -np.log10(df['pval_g1']+1e-6))

                ## output the statistic results
                variant_id = pval_g1 = beta_g1 = Z = pval_nominal = slope = z = pval_beta = is_eGene = np.nan
                if perm[perm['phenotype_id']==gene].shape[0]:
                    # gene exists in permutation file
                    variant_id = perm.loc[perm['phenotype_id']==gene, 'variant_id'].values[0]
                    if dflmm[dflmm['variant_id']==variant_id].shape[0]:
                        # top variants exist in lmm file
                        pval_g1 = dflmm.loc[dflmm['variant_id']==variant_id, 'pval_g1'].values[0]
                        beta_g1 = dflmm.loc[dflmm['variant_id']==variant_id, 'beta_g1'].values[0]
                        Z = dflmm.loc[dflmm['variant_id']==variant_id, 'Z'].values[0]
                        # results in permutation file
                        pval_nominal = perm.loc[perm['phenotype_id']==gene, 'pval_nominal'].values[0]
                        slope = perm.loc[perm['phenotype_id']==gene, 'slope'].values[0]
                        z = perm.loc[perm['phenotype_id']==gene, 'z'].values[0]
                        pval_beta = perm.loc[perm['phenotype_id']==gene, 'pval_beta'].values[0]
                        is_eGene = perm.loc[perm['phenotype_id']==gene, 'is_eGene'].values[0]

                ## output statistics
                print(f'{gene}: statistic writed.\n')
                f.write(f'{tissue}\t{gene}\t{rs}\t{rz}\t{rp}\t{cis_num}\t{variant_id}\t{pval_g1}\t{beta_g1}\t{Z}\t{pval_nominal}\t{slope}\t{z}\t{pval_beta}\t{is_eGene}\n')


if __name__ == '__main__':
    main()
