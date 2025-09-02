import vcfpy
import csv
import argparse

# 将基因型转换为数值
def genotype_to_numeric(genotype):
    if genotype in ['0/0', '0|0']:
        return 0  # Homozygous reference
    elif genotype in ['0/1', '1/0', '0|1', '1|0']:
        return 1  # Heterozygous
    elif genotype in ['1/1', '1|1']:
        return 2  # Homozygous variant
    else:
        return -1  # Unknown or missing genotype

# 读取VCF文件并转换为CSV
def vcf_to_encoded_csv(vcf_path, csv_path):
    # 读取VCF文件
    reader = vcfpy.Reader.from_path(vcf_path)

    # 创建CSV文件并写入标题
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # 写入标题行，包含所有样本名
        writer.writerow(['Position'] + reader.header.samples.names)

        # 遍历VCF文件中的每一个记录
        for record in reader:
            position = f"{record.CHROM}_{record.POS}"
            row = [position]
            for call in record.calls:
                genotype = call.data.get('GT', './.')
                numeric_genotype = genotype_to_numeric(genotype)
                row.append(numeric_genotype)
            writer.writerow(row)  # 写入编码后的基因型数据

# 主函数
def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="Convert a VCF file to an encoded CSV format.")
    parser.add_argument("vcf_path", help="Path to the input VCF file")
    parser.add_argument("csv_path", help="Path to the output CSV file")
    args = parser.parse_args()

    # 执行转换
    vcf_to_encoded_csv(args.vcf_path, args.csv_path)

if __name__ == "__main__":
    main()
