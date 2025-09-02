import os
import pandas as pd
import numpy as np

# 定义文件路径
base_dir = "/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split"
tissue_list_file = "/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/tissue40.list"
output_dir = "/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing"
output_file = os.path.join(output_dir, "PCGlnc_mat.txt")

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 从文件中读取组织列表
tissue_df = pd.read_csv(tissue_list_file, sep="\t", header=0)
tissue_list = tissue_df.iloc[:, 0].tolist()

combined_df = pd.DataFrame()

# 遍历组织列表
for tissue in tissue_list:
    file_path = os.path.join(base_dir, tissue, "phenotypes", f"{tissue}.expression.bed.gz")

    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"文件 {file_path} 不存在，跳过该组织")
        continue

    try:
        # 读取文件
        df = pd.read_csv(file_path, sep="\t", compression='gzip', header=0)

        # 打印调试信息
        print(f"处理文件: {file_path}")
        print(f"数据框形状: {df.shape}")

        # 提取第四列作为行名
        df.set_index(df.columns[3], inplace=True)

        # 计算第五列到最后一列的中位数
        median_values = df.iloc[:, 4:].median(axis=1)

        # 打印调试信息
        print(f"中位数计算结果: {median_values.head()}")

        # 创建一个新的数据框，包含组织名称和中位数
        tissue_df = pd.DataFrame(median_values, columns=[tissue])

        # 将结果合并到总的数据框中
        if combined_df.empty:
            combined_df = tissue_df
        else:
            combined_df = combined_df.join(tissue_df, how='outer')

    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {e}")

# 输出合并的矩阵到文件
combined_df = combined_df.sort_index(axis=1)
combined_df.index.name = ""
combined_df.to_csv(output_file, sep="\t")

print(f"合并的矩阵已保存到 {output_file}")
