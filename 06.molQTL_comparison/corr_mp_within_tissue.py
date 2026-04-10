# -*- coding: utf-8 -*-
"""
# @FileName      : corr_mp_within_tissue.py
# @Time          : 2024-07-19 19:03:35
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   :
"""

import click
import pandas as pd

@click.command()
@click.option('--infile1', help=".expressed.bed.gz file of mol phenotype1")
@click.option('--infile2', help=".expressed.bed.gz file of mol phenotype2")
@click.option('--outfile')
def main(infile1, infile2, outfile):
    """
    Calculate the pearson correlation between two different molecular phenotypes in a gene within a tissue.
    """
    tissue = infile1.strip().split('/')[0]
    df1 = pd.read_csv(infile1, sep="\t")
    df2 = pd.read_csv(infile2, sep="\t")

    if "gene_id" in df1.columns:
        df1 = df1.set_index(['gene_id', df1.columns[3]]).iloc[:, 3:]
    else: # phenotype id as gene id
        df1 = df1.set_index(df1.columns[3]).iloc[:, 3:].rename_axis(index='gene_id')
        
    if "gene_id" in df2.columns:
        df2 = df2.set_index(['gene_id', df2.columns[3]]).iloc[:, 3:]
    else: # phenotype id as gene id
        df2 = df2.set_index(df2.columns[3]).iloc[:, 3:].rename_axis(index='gene_id')

    common_column = df1.columns.intersection(df2.columns)
    df1 = df1[common_column].T
    df2 = df2[common_column].T

    corr = df1.corrwith(df2).dropna()
    df = corr.to_frame().reset_index()
    df.columns = [f"pairs", "Pearson correlation"]
    df['tissue'] = tissue
    df.to_csv(outfile, sep='\t', index=None)


if __name__ == '__main__':
    main()
