# -*- coding: utf-8 -*-
"""
# @FileName      : corr_mp_within_tissue.py
# @Time          : 2024-07-19 19:03:35
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   :
"""

import click
import numpy as np
import pandas as pd

@click.command()
@click.option('--infile', help=".expressed.bed.gz file of isoform, with colmuns of gene_id and phenotype_id")
@click.option('--outfile')
def main(infile, outfile):
    """
    Calculate the pearson correlation between different isoforms in a gene within a tissue.
    """
    tissue = infile.strip().split('/')[0]
    dfi = pd.read_csv(infile, sep="\t")
    dfi = dfi.set_index(['gene_id', dfi.columns[3]]).iloc[:, 3:]
    
    geneset = set(dfi.index.get_level_values('gene_id'))
    dfo = pd.DataFrame()
    for i, gene in enumerate(geneset, 1):
        df_ = dfi.loc[gene].T.corr()
        df_ = df_.where(np.triu(np.ones(df_.shape), 1).astype(np.bool_)) # get the upper triangle
        df_ = df_.reset_index().melt(id_vars='phenotype_id', value_vars=df_.columns, var_name='phenotype_id2', value_name='Pearson correlation').dropna().rename(columns={'phenotype_id':'phenotype_id1'})
        df_['gene_id'] = gene
        dfo = pd.concat([dfo, df_])
        print(f'{i}/{len(geneset)}')
        
    dfo['tissue'] = tissue
    dfo.to_csv(outfile, sep='\t', index=None)


if __name__ == '__main__':
    main()
