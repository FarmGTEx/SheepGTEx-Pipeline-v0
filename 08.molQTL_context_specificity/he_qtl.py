# -*- coding: utf-8 -*-
"""
# @FileName      : he_qtl.py
# @Time          : 2024-10-22 21:08:50
# @Author        : Mian Gong
# @Email         : gongmian2767@126.com
# @description   : Calculate he-QTL.
"""


import click
import pandas as pd
from scipy.stats import norm
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


@click.command()
@click.option('--susier_file', help="The input file of shared cs of susier results between two populations.")
@click.option('--mesusie_file', help="The input file of shared cs of mesusie results.")
@click.option('--outfile', help="The output file of he-eQTL.")
def main(susier_file, mesusie_file, outfile):
    # prepare data
    susier = pd.read_csv(susier_file, sep='\t')
    susier['PIP'] = (susier['PIP']+susier['PIP.1'])/2 # mean
    susier['type'] = 'susier'
    mesusie = pd.read_csv(mesusie_file, sep='\t')
    mesusie['type'] = 'mesusie'
    #filtered_mesusie = mesusie.groupby(['tissue', 'pheno_id', 'cs']).filter(lambda x: (x['EUR_CEA'] > 0.5).any())
    df = pd.concat([susier[['type','tissue','pheno_id','SNP','cs','PIP','MAF','BETA','SE','Z','N','MAF.1','BETA.1','SE.1','Z.1','N.1']],
                    mesusie[['type','tissue','pheno_id','SNP','cs','PIP','MAF','BETA','SE','Z','N','MAF.1','BETA.1','SE.1','Z.1','N.1']]])
    df = df[(df['MAF']>0.01)&(df['MAF.1']>0.01)]

    # clculate significance
    ngene = df.groupby('tissue')['pheno_id'].nunique().to_dict()
    df_sig = pd.DataFrame()
    for tissue in ngene:
        print(tissue)
        Q_value_threshold = 0.05/ngene[tissue]
        z = norm.ppf(1 - Q_value_threshold/2)
        df_ = df[df['tissue']==tissue].copy()
        df_['ngene'] = ngene[tissue]
        df_['Q_value_threshold'] = Q_value_threshold
        df_['CI1_lower'] = df_['BETA'] - z * df_['SE']
        df_['CI1_upper'] = df_['BETA'] + z * df_['SE']
        df_['CI2_lower'] = df_['BETA.1'] - z * df_['SE.1']
        df_['CI2_upper'] = df_['BETA.1'] + z * df_['SE.1']
        df_['sig'] = (df_['CI1_lower'] > df_['CI2_upper']) | (df_['CI2_lower'] > df_['CI1_upper'])
        df_sig = pd.concat([df_sig, df_])

    # save results
    df_sig.to_csv(outfile, index=None, sep="\t")
    

if __name__ == '__main__':
    main()
