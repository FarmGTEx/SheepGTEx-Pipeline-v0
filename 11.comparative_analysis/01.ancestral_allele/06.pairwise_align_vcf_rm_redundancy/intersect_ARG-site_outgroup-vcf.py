#baifengting
#20230331

import gzip
import cyvcf2
import click

@click.command()
@click.option('--vcf_arg', required=True, help='用于构建ARG的vcf文件，用于过滤outgroup vcf文件中的位点')
@click.option('--raw_vcf', required=True, help='原始outgroup vcf')
@click.option('--out_vcf')
def main(vcf_arg, raw_vcf, out_vcf):
    st_arg_site = set()
    vcf_ARG = cyvcf2.VCF(vcf_arg)
    for variant in vcf_ARG:
        st_arg_site.add(variant.POS)
    with gzip.open(raw_vcf, 'r') as f_in, gzip.open(out_vcf, 'w') as f_out:
        for line in f_in:
            line_decode = line.decode().strip()
            if line_decode.startswith('#'):
                f_out.write(line)
            else:
                ls_line = line_decode.split()
                if int(ls_line[1]) in st_arg_site:
                    f_out.write(line)

if __name__ == '__main__':
    main()