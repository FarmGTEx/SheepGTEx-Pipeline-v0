#baifengting
#20230508

import gzip
import numpy as np
import cyvcf2
import click

def get_ref_info(vcf_ref):
    #return dic_ref_gt = {pos1:gt1, pos2:...}, st_ref_pos
    base_index_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    dic_ref_gt = {}
    st_ref_pos = set()
    vcf = cyvcf2.VCF(vcf_ref)
    total_alleles = len(vcf.samples) * 2
    for variant in vcf:
        ay_gt = np.zeros(4, dtype=int)
        pos = variant.POS
        ref_a = variant.REF
        alt_a = variant.ALT[0]
        num_alt = variant.num_het + variant.num_hom_alt * 2
        num_ref = total_alleles - num_alt
        ay_gt[base_index_map[ref_a]] = num_ref
        ay_gt[base_index_map[alt_a]] = num_alt
        gt = '{:},{:},{:},{:}'.format(ay_gt[0], ay_gt[1], ay_gt[2], ay_gt[3])
        dic_ref_gt[pos] = gt
        st_ref_pos.add(pos)
    return dic_ref_gt, st_ref_pos

def get_est_gt(vcf_outgroup, index, st_ref_pos, allele_input_map, dic_input_gt, dic_ref_gt):
    st_pos = set(dic_input_gt.keys())
    with gzip.open(vcf_outgroup, 'r') as f_out:
        for line in f_out:
            line = line.decode().strip()
            if line.startswith('#'):
                pass
            else:
                line = line.split()
                if 'PASS' in [line[5], line[6]]:
                    pos = int(line[1])
                    assert pos in st_ref_pos #输入的outgroup file应只包含用于构建ARG的位点
                    ref_a = line[3]
                    if line[4] == '.' :
                        alt_a = ref_a
                    else:
                        alt_a = line[4]
                        #if pos not in set(dic_input_gt.keys()):
                    if pos not in st_pos:
                        ls_gt = [dic_ref_gt[pos]] + ['0,0,0,0'] * 3
                        ls_gt[index] = allele_input_map[alt_a]
                        dic_input_gt[pos] = ls_gt
                        st_pos.add(pos)
                    else:
                        assert dic_input_gt[pos][index] == '0,0,0,0' # 检验此位点没有重复基因型
                        dic_input_gt[pos][index] = allele_input_map[alt_a]
                else:
                    pass

@click.command()
@click.option('--vcf_ref', required=True, help='focal species的代表性群体组成的vcf，比如千牛中每个品种各选1个高深度的个体')
@click.option('--vcf_outgroup1', required=True)
@click.option('--vcf_outgroup2', required=True)
@click.option('--vcf_outgroup3', required=True)
@click.option('--outfile', required=True)
def main(vcf_ref, vcf_outgroup1, vcf_outgroup2, vcf_outgroup3, outfile):
    """
    此脚本中，输入的vcf_ref 需要与用于构建ARG的vcf文件位点一致，且外群vcf的位点仅保留了同时存在于vcf_ref中的位点，vcf_outgroup = vcf_outgroup & vcf_ref.\n
    输出文件为esf-sfs输入文件outfile以及每个位点物理坐标的outfile.site
    """
    allele_input_map = {
    'A': '1,0,0,0', 
    'C': '0,1,0,0', 
    'G': '0,0,1,0', 
    'T': '0,0,0,1', 
    }

    dic_ref_info, st_ref_pos = get_ref_info(vcf_ref)
    dic_input_gt = {}
    for index, vcf_outgroup in enumerate([vcf_outgroup1, vcf_outgroup2, vcf_outgroup3], start=1):
        get_est_gt(vcf_outgroup, index, st_ref_pos, allele_input_map, dic_input_gt, dic_ref_info)
    
    with open(outfile, 'w') as f_out, open(outfile + '.site', 'w') as f_site:
        for pos, ls_gt in dic_input_gt.items():
            f_site.write(str(pos) + '\n')
            f_out.write(' '.join(ls_gt) + '\n')

if __name__ == '__main__':
    main()