"""Assemble data into an RNA phenotype BED file"""
"""Reused from https://github.com/PejLab/Pantry/blob/main/Project/scripts/assemble_bed.py"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import re

def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations from GTF file

    Returns TSS as the first four columns of the BED format, meaning the
    coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    gene_tss_positions = []

    with open(ref_anno, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            chrom = columns[0]
            feature_type = columns[2]
            if feature_type != 'gene':
                continue

            attributes = columns[8]
            if 'gene_id' not in attributes:
                continue
            gene_id = re.findall(r'gene_id "(.+?)"', attributes)[0]

            if gene_id is None:
                continue

            strand = columns[6]
            if strand == '+':
                tss_position = int(columns[3])
            elif strand == '-':
                tss_position = int(columns[4])
            else:
                continue

            gene_tss_positions.append([chrom, tss_position-1, tss_position, gene_id])

    # Name columns for tensorQTL:
    anno = pd.DataFrame(gene_tss_positions, columns=['#chr', 'start', 'end', 'gene_id'])
    return anno

def load_featureCounts(sample_ids: list, counts_dir: Path, feature: str, min_count: int = 10) -> pd.DataFrame:
    """Assemble featureCounts outputs into a table"""
    counts = []
    for i, sample in enumerate(sample_ids):
        fname_ex = counts_dir / f'{sample}.{feature}.tsv'
        d = pd.read_csv(fname_ex, sep='\t', index_col='Geneid', skiprows=1)
        d = d.iloc[:, 5]
        d.name = sample
        d[d < min_count] = np.nan
        # Store in lists and concat all at once to avoid 'PerformanceWarning: DataFrame is highly fragmented' warning
        counts.append(d)
    return pd.concat(counts, axis=1)

def assemble_stability(sample_ids: list, stab_dir: Path, ref_anno: Path, bed: Path):
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    print("Load exon counts...")
    exon = load_featureCounts(sample_ids, stab_dir, 'consExons')
    print("Load intron counts...")
    intron = load_featureCounts(sample_ids, stab_dir, 'introns')
    genes = exon.index[np.isin(exon.index, intron.index)]
    # genes = set(exon.index).intersection(intron.index)
    assert exon.loc[genes, :].index.equals(intron.loc[genes, :].index)
    assert exon.columns.equals(intron.columns)
    print("Calculate stability and filtering genes...")
    df = exon.loc[genes, :] / intron.loc[genes, :]
    df = df[df.isnull().mean(axis=1) <= 0.5] # remove genes with more than half missing value
    print("Load tss...")
    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'gene_id': 'phenotype_id'})
    df.index = df.index.rename('phenotype_id')
    df = anno.merge(df.reset_index(), on='phenotype_id', how='inner')
    print("Save unnormalized results...")
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

parser = argparse.ArgumentParser(description='Assemble data into an RNA phenotype BED file')
parser.add_argument('--input-dir', type=Path, required=True, help='Directory containing input files, for phenotype groups with per-sample input files')
parser.add_argument('--samples', type=Path, required=True, help='File listing sample IDs, for phenotype groups with per-sample input files')
parser.add_argument('--ref_anno', type=Path, required=True, help='Reference annotation file')
parser.add_argument('--output', type=Path, required=True, help='Output file ("*.bed")')
args = parser.parse_args()

samples = pd.read_csv(args.samples, sep='\t', header=None, dtype=str)[0].tolist()
assemble_stability(samples, args.input_dir, args.ref_anno, args.output)
