from __future__ import print_function
#import datatable as dt
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import gzip
import contextlib
from datetime import datetime
import tempfile
import shutil
import glob
from sklearn.decomposition import PCA

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def get_ann_bed(annotation_bed, feature='gene'):
    """
    Parse genes from GTF, create placeholder DataFrame for BED output
    """
    chrom = []
    start = []
    end = []
    gene_id = []
    with open(annotation_bed, 'r') as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            chrom.append(row[0])

            # TSS: gene start (0-based coordinates for BED)
            if row[3]=='+':
                start.append(np.int64(row[1])-1)
                end.append(np.int64(row[1]))
            elif row[3]=='-':
                start.append(np.int64(row[2])-1)  # last base of gene
                end.append(np.int64(row[2]))
            else:
                raise ValueError('Strand not specified.')

            gene_id.append(row[4])

    bed_df = pd.DataFrame(
        data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id},
        columns=['chr', 'start', 'end', 'gene_id'],
        index=gene_id)
    return bed_df


def write_bed(bed_df, output_name):
    """Write DataFrame to BED"""
    bgzip = subprocess.Popen('bgzip -c > '+output_name,
        stdin=subprocess.PIPE, shell=True)
    bed_df.to_csv(bgzip.stdin, sep='\t', index=False)
    stdout, stderr = bgzip.communicate()
    subprocess.check_call('tabix -f '+output_name, shell=True)


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run leafcutter clustering, prepare for FastQTL')
    parser.add_argument('sample_participant_lookup', help='Lookup table linking samples to participants')
    #parser.add_argument('exons', help='Exon definitions file, with columns: chr, start, end, strand, gene_id, gene_name')
    parser.add_argument('genes_bed', help='Collapsed gene annotation in bed format')
    parser.add_argument('inprefix', help='Prefix for input files (sample set ID)')
    parser.add_argument('outprefix', help='Prefix for output files (sample set ID)')
    parser.add_argument('--num_pcs', default=10, type=int, help='Number of principal components to calculate')
    parser.add_argument('--leafcutter_dir', default='/opt/leafcutter',
                        help="leafcutter directory, containing 'clustering' directory")
    parser.add_argument('-i', '--input_dir', default='.', help='Input directory')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('-c', '--chunksize', default=100000, type=int, help='chunk size for the input file')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] leafcutter clustering')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup, sep='\t', index_col=0, header=None, dtype=str).squeeze('columns')

    filtered_counts_file = os.path.join(args.output_dir, args.outprefix+'_perind.counts.filtered.gz')
    all_counts_df = pd.read_csv(os.path.join(args.input_dir, args.inprefix+'_perind.counts.gz'), sep='\s+', index_col=0, chunksize=args.chunksize)
    for i, counts_df in enumerate(all_counts_df):
        # change columns to sample IDs
        col_dict = {i:i.split('_Aligned.sortedByCoord.out.filtered.bam')[0] for i in counts_df.columns}
        counts_df.rename(columns=col_dict, inplace=True)
        assert sample_participant_lookup_s.index.isin(counts_df.columns).all()
        # filter samples
        counts_df = counts_df[sample_participant_lookup_s.index.to_list()]
    
        print(f'  * filtering counts in chunk {i}')
        calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1]>0 else 0
        frac_df = counts_df.applymap(lambda x: calculate_frac([int(i) for i in x.split('/')]))
        pct_zero = (frac_df==0).sum(axis=1) / frac_df.shape[1]  # for zero counts, frac is zero
        n_unique = frac_df.apply(lambda x: len(x.unique()), axis=1)
        zscore_df = ((frac_df.T-frac_df.mean(1)) / frac_df.std(1)).T

        # filter out introns with low counts or low complexity
        n = np.floor(frac_df.shape[1]*0.1)
        if n<10:
            n = 10
        mask = (pct_zero <= 0.5) & (n_unique >= n)
        # additional filter for low complexity
        ns = zscore_df.shape[1]
        mask2 = ((zscore_df.abs()<0.25).sum(1)>=ns-3) & ((zscore_df.abs()>6).sum(1)<=3)
        if np.any(mask & mask2):
            print('    ** dropping {} introns with low variation'.format(np.sum(mask & mask2)))
        mask = mask & ~mask2

        filtered_counts_df = counts_df.loc[mask].copy()
        cluster_ids = np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))
        filtered_cluster_ids = np.unique(filtered_counts_df.index.map(lambda x: x.split(':')[-1]))
        print('    ** dropping {} introns with counts in fewer than 50% of samples\n'
                '       {}/{} introns remain ({}/{} clusters)'.format(
                    counts_df.shape[0]-filtered_counts_df.shape[0], filtered_counts_df.shape[0],
                    counts_df.shape[0], len(filtered_cluster_ids), len(cluster_ids)))
        print(f'  * writing counts in chunk {i}')
        if i == 0:
            filtered_counts_df.to_csv(filtered_counts_file, sep=' ', mode='w')
        else:
            filtered_counts_df.to_csv(filtered_counts_file, header=None, sep=' ', mode='a')

    print('  * preparing phenotype table')
    subprocess.check_call(
        '/storage/public/home/2022060207/anaconda3/envs/leafcutter/bin/python2 '+os.path.join(args.leafcutter_dir, 'scripts', 'prepare_phenotype_table.py') \
        +' '+filtered_counts_file \
        +' -p '+str(args.num_pcs), shell=True)

    print('  * concatenating chromosome-level BED files')
    bed_files = sorted(glob.glob(os.path.join(args.output_dir, args.outprefix+'_perind.counts.filtered.gz.qqnorm_*')))
    bed_df = []
    for f in bed_files:
        # only remain autosomes
        if ('qqnorm_NW_' not in f) and ('qqnorm_chrX' not in f) and ('qqnorm_chrY' not in f) and ('qqnorm_chrMT' not in f):
            bed_df.append(pd.read_csv(f, sep='\t', dtype=str))
    bed_df = pd.concat(bed_df, axis=0)

    print('    ** sorting')
    # leafcutter doesn't produce output for chrX --> numeric sort
    bed_df['chr_ix'] = bed_df['#Chr'].str.replace('chr','').str.replace('X','100').astype(np.int32)
    for c in ['start', 'end']:
        bed_df[c] = bed_df[c].astype(np.int32)
    bed_df = bed_df.sort_values(['chr_ix', 'start', 'end'],ascending = (True, True, True))
    bed_df.drop('chr_ix', axis=1, inplace=True)

    print('    ** writing BED')
    bed_file = os.path.join(args.output_dir, args.outprefix+'.perind.counts.filtered.qqnorm.bed.gz')
    bgzip = subprocess.Popen('bgzip -c -f > '+bed_file, stdin=subprocess.PIPE, shell=True, universal_newlines=True)
    stdout, stderr = bgzip.communicate(bed_df.to_csv(sep='\t', index=False))
    print('    ** indexing')
    subprocess.check_call('tabix -f '+bed_file, shell=True)

    print('  * converting cluster coordinates to gene coordinates')
    tss_df = get_ann_bed(args.genes_bed)
    cluster2gene_dict = pd.read_csv(os.path.join(args.input_dir, args.inprefix+'.leafcutter.clusters_to_genes.txt'), sep='\t', index_col=0, low_memory=False).squeeze('columns').to_dict()

    print('    ** assigning introns to gene mapping(s)')
    n = 0
    gene_bed_df = []
    group_s = {}
    for _,r in bed_df.iterrows():
        s = r['ID'].split(':')
        #cluster_id = s[0]+':'+s[-1]
        cluster_id = r['ID']
        if cluster_id in cluster2gene_dict:
            gene_ids = cluster2gene_dict[cluster_id].split(',')
            for g in gene_ids:
                gi = r['ID']+':'+g
                gene_bed_df.append(tss_df.loc[g, ['chr', 'start', 'end']].tolist() + [gi] + r.iloc[4:].tolist() + [g])
                group_s[gi] = g
        else:
            n += 1
    if n>0:
        print('    ** discarded {} introns without gene mapping'.format(n))
    gene_bed_df = pd.DataFrame(gene_bed_df, columns=bed_df.columns.tolist()+['gene_id'])

    print('    ** sorting for gene id')
    gene_bed_df['chr_ix'] = gene_bed_df['#Chr'].str.replace('chr','').str.replace('X','100').astype(np.int32)
    for c in ['start', 'end']:
        gene_bed_df[c] = gene_bed_df[c].astype(np.int32)
    gene_bed_df = gene_bed_df.sort_values(['chr_ix', 'start', 'end', 'gene_id'],ascending = (True, True, True, True))
    gene_bed_df.drop(['chr_ix', 'gene_id'], axis=1, inplace=True)
    
    print('  * writing tensorQTL inputs')
    # change sample IDs to participant IDs
    gene_bed_df.rename(columns=sample_participant_lookup_s, inplace=True)
    #write_bed(gene_bed_df, os.path.join(args.output_dir, args.prefix+'.leafcutter.bed.gz'))
    gene_bed_file=os.path.join(args.output_dir, args.outprefix+'.leafcutter.sorted.bed')
    gene_bed_df.to_csv(gene_bed_file, sep='\t', index=False)
    subprocess.check_call('bgzip -f '+gene_bed_file, shell=True)
    subprocess.check_call('tabix -p bed -f '+gene_bed_file+'.gz', shell=True)
    pd.Series(group_s, dtype='object').sort_values().to_csv(os.path.join(args.output_dir, args.outprefix+'.leafcutter.phenotype_groups.txt'), sep='\t', header=False)

    print('  * calculating PCs')
    pca = PCA(n_components=args.num_pcs)
    pca.fit(bed_df[bed_df.columns[4:]])
    pc_df = pd.DataFrame(pca.components_, index=['PC{}'.format(i) for i in range(1,11)],
        columns=['-'.join(i.split('-')[:2]) for i in bed_df.columns[4:]])
    pc_df.index.name = 'ID'
    pc_df.to_csv(args.output_dir+'/'+args.outprefix+'.leafcutter.PCs.txt', sep='\t')

print('done')
