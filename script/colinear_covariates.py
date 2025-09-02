#!/usr/bin/env python3
# Author: Mian Gong
import pandas as pd
import numpy as np
import click


def category(covardf):
    for Covariates in covardf.columns:
        covarset = set(covardf[Covariates])
        if len(covarset) == 1: # duplicated data
            print(f"Drop constant covariates: {Covariates}")
            del covardf[Covariates]
            continue
        if covardf[Covariates].dtype == 'O': # category data
            if len(covarset) >= 2:
                removed_element = covarset.pop()  # remove an element to make a non-singular matrix
                print("Removed element:", removed_element)
                for v in covarset:
                    covardf[v] = 0
                    covardf.loc[covardf[Covariates]==v, v] = 1
                print(f"Split multi-category covariate into 2-category covariates: {Covariates}")
                del covardf[Covariates]
    covardf = covardf.astype(float)
    return covardf

def drop_colinear(combined_df):
    # Method from https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/combine_covariates.py
    C = combined_df.astype(np.float64)
    Q,R = np.linalg.qr(C-np.mean(C, axis=0))
    colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
    if np.any(colinear_ix):
        print('Colinear covariates detected:')
        for i in C.columns[colinear_ix]:
            print("  * dropped '{}'".format(i))
        combined_df = combined_df.loc[:,~colinear_ix]

    # drop colinear covariates with R2 >= 0.9
    dropped = []
    for i in range(combined_df.shape[1] - 1):
        col1 = combined_df.iloc[:, i]
        for j in range(i + 1, combined_df.shape[1]):
            col2 = combined_df.iloc[:, j]
            r2 = col1.corr(col2) ** 2
            if (col1.name not in dropped) and (col1.name not in dropped) and (r2 >= 0.9):
                print(f"R2 between columns {col1.name} and {col2.name}: {r2}. Removing column {col2.name}")
                dropped.append(col2.name)
    combined_df = combined_df.drop(columns=dropped)
    return combined_df


@click.command()
@click.option('--outfile')
@click.argument('covarfiles', nargs=-1)
def main(outfile, covarfiles):
    '''
    Identify and drop colinear covariates. The first column of each input file is sample name.
    '''
    for i, covarfile in enumerate(covarfiles):
        if i == 0:
            combined_df = category(pd.read_csv(covarfile, sep='\t', index_col=0))
        else:
            additional_df = category(pd.read_csv(covarfile, sep='\t', index_col=0))
            combined_df = pd.merge(combined_df, additional_df, left_index=True, right_index=True)

    # identify and drop colinear covariates
    combined_df = drop_colinear(combined_df)
    combined_df.to_csv(outfile, sep='\t')

if __name__ == '__main__':
    main()
    
