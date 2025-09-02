from __future__ import print_function
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


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run leafcutter clustering, prepare for FastQTL')
    parser.add_argument('junc_files_list', help='File with paths to ${sample_id}.junc files')
    parser.add_argument('exons', help='Exon definitions file, with columns: chr, start, end, strand, gene_id, gene_name')
    #parser.add_argument('genes_bed', help='Collapsed gene annotation in bed format')
    parser.add_argument('prefix', help='Prefix for output files (sample set ID)')
    parser.add_argument('--min_clu_reads', default='50', type=str, help='Minimum number of reads supporting each cluster')
    parser.add_argument('--min_clu_ratio', default='0.001', type=str, help='Minimum fraction of reads in a cluster that support a junction')
    parser.add_argument('--max_intron_len', default='500000', type=str, help='Maximum intron length')
    parser.add_argument('--leafcutter_dir', default='/opt/leafcutter',
                        help="leafcutter directory, containing 'clustering' directory")
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] leafcutter clustering')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)


    print('  * running leafcutter clustering')
    # generates ${prefix}_perind_numers.counts.gz and ${prefix}_perind.counts.gz
    subprocess.check_call(
        '/storage/public/home/2022060207/anaconda3/envs/leafcutter/bin/python2 '+os.path.join(args.leafcutter_dir, 'clustering', 'leafcutter_cluster.py' \
            +' --juncfiles '+args.junc_files_list \
            +' --outprefix '+args.prefix \
            +' --rundir '+args.output_dir \
            +' --minclureads '+args.min_clu_reads \
            +' --mincluratio '+args.min_clu_ratio \
            +' -s True --maxintronlen '+args.max_intron_len), shell=True)

    print('  * compressing outputs')
    subprocess.check_call('gzip -f '+os.path.join(args.output_dir,args.prefix+'_pooled'), shell=True)
    subprocess.check_call('gzip -f '+os.path.join(args.output_dir,args.prefix+'_refined'), shell=True)

    print('  * mapping clusters to genes')
    subprocess.check_call(
        'Rscript' \
            +' map_clusters_to_genes.R' \
            +' '+os.path.join(args.output_dir, args.prefix+'_perind.counts.gz') \
            +' '+args.exons \
            +' '+args.output_dir+'/'+args.prefix + '.leafcutter.clusters_to_genes.txt', shell=True)
print('done')
