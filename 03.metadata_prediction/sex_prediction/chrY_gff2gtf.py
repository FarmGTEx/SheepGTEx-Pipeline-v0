# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 18:20:46 2023

@author: Gongmian

E-mail: gongmian2767@126.com
"""

import click
import re
from collections import defaultdict

@click.command()
@click.option('--ingff')
@click.option('--outgtf')
def main(ingff, outgtf):
    with open(ingff, 'r') as fi, open(outgtf, 'w') as fo:
        tid2gid = dict() ; exonum = defaultdict(int)
        for line in fi:
            tline = line.strip().split()
            tline[0] = 'chrY'
            #if len(tline)>1 and tline[2] == 'exon':
            #    tid = re.findall("Parent=\S*", tline[8])[0].split('=')[1].replace(';','')
            #    exonum[tid] += 1
            #    fo.write('\t'.join(tline[0:8]))
            #    fo.write(f'\tgene_id "{tid2gid[tid]}"; transcript_id "{tid}"; gene "{tid2gid[tid]}"; transcript_biotype "mRNA"; exon_number "{exonum[tid]}";\n')
            if len(tline)>1 and tline[2] == 'gene':
                tline[2] = tline[2].replace('gene','transcript')
                tid = re.findall("ID=\S*?;", tline[8])[0].split('=')[1].replace(';','')
                gid = re.findall("Target=\S*", tline[8])[0].split('=')[1].replace(';','')
                tid2gid[tid] = gid
                fo.write('\t'.join(tline[0:8]))
                fo.write(f'\tgene_id "{tid2gid[tid]}"; transcript_id "{tid}"; gbkey "mRNA"; gene "{tid2gid[tid]}"; transcript_biotype "mRNA";\n')
            elif len(tline)>1 and tline[2] == 'mRNA':
                tline[2] = tline[2].replace('mRNA','exon')
                tid = re.findall("ID=\S*?;", tline[8])[0].split('=')[1].replace(';','')
                gid = re.findall("Target=\S*", tline[8])[0].split('=')[1].replace(';','')
                exonum[tid] += 1
                fo.write('\t'.join(tline[0:8]))
                fo.write(f'\tgene_id "{tid2gid[tid]}"; transcript_id "{tid}"; gene "{tid2gid[tid]}"; transcript_biotype "mRNA"; exon_number "{exonum[tid]}";\n')
                
if __name__ == '__main__':
    main()

