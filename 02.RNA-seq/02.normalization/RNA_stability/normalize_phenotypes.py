"""Normalize phenotypes in a BED file for direct use in QTL mapping"""
"""Reused from https://github.com/PejLab/Pantry/blob/main/Project/scripts/normalize_phenotypes.py"""

import argparse, os
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.impute import KNNImputer


def impute_knn(X):
    """Impute missing values using KNN imputation"""
    imputer = KNNImputer(n_neighbors=10, weights='uniform')
    return imputer.fit_transform(X)


def normalize_quantiles(df):
    """
    Quantile normalization to the average empirical distribution
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")
    Author: Francois Aguet

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    M = df.values.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=int)  # replaced deprecated np.int with int -D. Munro
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    Author: Francois Aguet
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


parser = argparse.ArgumentParser(description="Normalize phenotypes in a BED file for direct use in QTL mapping")
parser.add_argument("-i", "--input", help="Input BED file", required=True)
parser.add_argument("-s", "--samples", help="TSV linking sample IDs (columns in expression matrices) to participant IDs (VCF IDs)", required=False)
parser.add_argument("-o", "--output", help="Output BED file", required=True)
parser.add_argument("-c", "--chrs", help='File listing chromosomes to include (default: chr1-22 + chrX)', required=False)
# parser.add_argument("--upper-quartile", help="Prior to quantile transformation, scale values per sample such that the upper quartile value is the same across samples", action="store_true", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", dtype={"#chr": str, "start": int, "end": int, "phenotype_id": str})

# Stop if the raw table is empty:
assert not df.empty, f"{args.input} is empty"

if args.samples is not None:
    sample_to_participant_s = pd.read_csv(args.samples, sep='\t', index_col=0,
                                          header=None, dtype=str).squeeze('columns')
    missing_ids = ~sample_to_participant_s.index.isin(df.columns[4:])
    if missing_ids.any():
        raise ValueError(f"Some samples missing from expression file.")
    # filter samples
    df = df[list(df.columns[:4]) + sample_to_participant_s.index.to_list()]
    # change sample IDs to participant IDs'
    df.rename(columns=sample_to_participant_s, inplace=True)

if args.chrs is not None:
    with open(args.chrs) as f:
        chrs = f.read().strip().split('\n')
    df = df[df['#chr'].isin(chrs)]

data = df.drop(df.columns[:4], axis=1)

# Impute missing values, which are not allowed in tensorQTL
if data.isnull().values.any():
    print(f"{args.output}: {data.isnull().values.sum()} missing values")
    # remove rows with >50% missing values:
    mostly_missing = data.apply(lambda x: np.mean(x.isnull()) >= 0.5, axis=1)
    df = df[~mostly_missing]
    data = data[~mostly_missing]
    print(f"{args.output}: Removed {mostly_missing.sum()} rows with at least 50% missing values")
    print(f"{args.output}: Imputing {data.isnull().values.sum()} remaining missing values using K-nearest neighbors")
    data[:] = impute_knn(data.T).T

# Remove rows with more than 50% zeros, as recommended for tensorQTL:
mostly_zeros = data.apply(lambda x: np.mean(x == 0) > 0.5, axis=1)
df = df[~mostly_zeros]
data = data[~mostly_zeros]
# if mostly_zeros.sum() > 0:
print(f"{args.output}: Removed {mostly_zeros.sum()} rows with more than 50% zeros")

# Remove rows with no variation:
no_var = data.apply(lambda x: x.nunique() == 1, axis=1)
df = df[~no_var]
data = data[~no_var]
# if no_var.sum() > 0:
print(f"{args.output}: Removed {no_var.sum()} additional rows with no variation")

# Stop if the processed table is empty:
assert not df.empty, f"After filtering {args.input}, no phenotypes remain"

iqn = normalize_quantiles(data)
iqn = inverse_normal_transform(iqn)
# Replace data values in df with transformed values:
df.loc[:, data.columns] = iqn

print('    ** sorting')
df['chr_ix'] = df['#chr'].str.replace('chr','').str.replace('X','100').astype(np.int32)
for c in ['start', 'end']:
    df[c] = df[c].astype(np.int32)
df = df.sort_values(['chr_ix', 'start', 'end'],ascending = (True, True, True))
df.drop('chr_ix', axis=1, inplace=True)

df.to_csv(args.output, sep="\t", index=False, float_format="%g")
