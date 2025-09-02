# -*- coding: utf-8 -*-
"""
# @FileName      : tissue_specific_effect_qtl.py
# @Time          : 2024-10-22 21:08:50
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   : Calculate effect size specific-QTL between shared tissues.
"""

import click
import pyreadr
import numpy as np
import pandas as pd
from scipy.stats import norm
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
    

@click.command()
@click.option('--inlfsr', help="The input RDS file of lfsr score output by mashr.")
@click.option('--infile', help="The input file of joint pairwise tissue lead variants.")
@click.option('--outfile', help="The output file of he-eQTL.")
def main(inlfsr, infile, outfile):
    # prepare LFSR input
    lfsr = pyreadr.read_r(inlfsr)[None]
    lfsr_sig = (lfsr<0.05).reset_index()
    lfsr_sig['variant_id'] = lfsr_sig['index'].str.split(',', expand=True)[0]
    lfsr_sig['phenotype_id'] = lfsr_sig['index'].str.split(',', expand=True)[1]

    # prepare data
    df = pd.read_csv(infile, sep="\t")
    df = df[(df['is_eGene1']==True)&(df['is_eGene2']==True)]
    if df.empty:
        warnings.warn("No eGene exists in both tissues!", UserWarning)
    else:
        df['r2_min'] = df['r2'].apply(lambda x: min([float(i) for i in x.split(',') if i.strip()] ) if pd.notna(x) else np.nan)
        df.loc[df['variant_id1']==df['variant_id2'], 'r2_min'] = 1

        # pairwise shared significance
        tis1 = df['tissue1'].values[0]
        tis2 = df['tissue2'].values[0]
        ps_sig = lfsr_sig[['phenotype_id', 'variant_id', tis1, tis2]]
        ps_sig['shared'] = ps_sig[[tis1, tis2]].sum(axis=1)==2

        # merge
        df = pd.merge(df, ps_sig[['phenotype_id', 'variant_id', 'shared']], left_on=['phenotype_id', 'variant_id1'],
                      right_on=['phenotype_id', 'variant_id']).rename(columns={'shared':'shared1'})
        df = pd.merge(df, ps_sig[['phenotype_id', 'variant_id', 'shared']], left_on=['phenotype_id', 'variant_id2'], 
                      right_on=['phenotype_id', 'variant_id']).rename(columns={'shared':'shared2'})
    
        # clculate significance of effect size specific eQTLs
        #df_filtered = df[(df['r2_min']>=0.5)&(df['PP.H4.abf']>=0.5)&df['shared1']&df['shared2']]
        #df_filtered = df[(df['r2_min']>=0.8)&df['shared1']&df['shared2']]
        df_filtered = df
        #genes = df_filtered['phenotype_id'].values
        #Q_value_threshold = 0.05/(len(genes))
        Q_value_threshold = 0.05
        z = norm.ppf(1 - Q_value_threshold/2)

        ## get lead effect size confidence intervals per gene
        df_filtered['ci1_lower'] = df_filtered['slope1'] - z*df_filtered['slope_se1']
        df_filtered['ci1_higher'] = df_filtered['slope1'] + z*df_filtered['slope_se1']
        df_filtered['ci2_lower'] = df_filtered['slope2'] - z*df_filtered['slope_se2']
        df_filtered['ci2_higher'] = df_filtered['slope2'] + z*df_filtered['slope_se2']

        ## calculate significance
        df_filtered['type'] = 'Homogeneous'
        significant_effect = (df_filtered['ci1_lower']>df_filtered['ci2_higher'])|(df_filtered['ci2_lower']>df_filtered['ci1_higher'])
        df_filtered.loc[significant_effect, 'type'] = 'Heterogeneous'

        # get genes with opposite lead effect size between tissues
        opposite_effect = ((df_filtered['ci1_higher']<0)&(df_filtered['ci2_lower']>0))|((df_filtered['ci2_higher']<0)&(df_filtered['ci1_lower']>0))
        df_filtered.loc[opposite_effect, 'type'] = 'Opposite'

        # save results
        df_filtered.to_csv(outfile, index=None, sep="\t")

if __name__ == '__main__':
    main()
