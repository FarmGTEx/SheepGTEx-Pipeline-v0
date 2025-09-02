import pandas as pd
from scipy.stats import fisher_exact
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import argparse


# 定义计算OR值和95%置信区间的函数
def calculate_or_and_ci(a, b, c, d):
    # 检查零值
    if a == 0 or b == 0 or c == 0 or d == 0:
        or_value = float('inf') if (a * d > b * c) else 0
        ci_low, ci_high = (0, float('inf'))
        p_value = 1.0  # 如果不能计算，则P值设为1

    else:
        # 计算OR值
        or_value = (a * d) / (b * c)
    
        # 使用scipy计算P值
        table = [[a, b], [c, d]]
        _, p_value = fisher_exact(table, alternative='two-sided')
    
        # 使用statsmodels计算置信区间
        #table2x2 = sm.stats.Table2x2(np.array(table))
        #ci_low, ci_high = table2x2.oddsratio_confint()

        #计算标准误差和置信区间
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        ci_low= np.exp(np.log(or_value)) - np.exp( 1.96 * se)
        ci_high = np.exp(np.log(or_value)) + np.exp( 1.96 * se)
    
    return or_value, ci_low, ci_high, p_value

# 定义计算fold enrichment值和95%置信区间的函数
def calculate_fold_enrichment(target_count, target_total, background_count, background_total):
    # 检查零值
    if target_count == 0 or target_total == 0 or background_count == 0 or background_total == 0:
        fe_value = float('inf') if (target_count * background_total > background_count * target_total) else 0
        ci_low, ci_high = (0, float('inf'))
        p_value = 1.0  # 如果不能计算，则P值设为1
    
    else:
        # 计算目标组中某一项的比例
        target_ratio = target_count / target_total
        # 计算背景组中相同项的比例
        background_ratio = background_count / background_total
        # 计算 Fold Enrichment
        fe_value = target_ratio / background_ratio
    
        # 使用scipy计算P值
        table = [[target_count, target_total - target_count],
             [background_count, background_total - background_count]]
        _, p_value = fisher_exact(table, alternative='two-sided')
    
        # 使用statsmodels计算置信区间
        #table2x2 = sm.stats.Table2x2(np.array(table))
        #ci_low, ci_high = table2x2.oddsratio_confint()

        # 计算标准误差和置信区间
        se = np.sqrt(1/target_count + 1/target_total + 1/background_count + 1/background_total)
        ci_low = np.exp(np.log(fe_value)) - np.exp( 1.96 * se)
        ci_high = np.exp(np.log(fe_value)) + np.exp( 1.96 * se)
    
    return fe_value, ci_low, ci_high, p_value


# 文件目录设置
parser = argparse.ArgumentParser(description='calculate the OR of QTL and finemapping sites.')
parser.add_argument('qtl_dir', help='directory of files of xQTLs and finemapping sites')
parser.add_argument('annotation_dir', help='directory of files of functional annotation sites.')
parser.add_argument('site_file', help='file of all sites.')
parser.add_argument('out_dir', help='directory to save results.')
args = parser.parse_args()

# 获取目录中的文件列表
qtl_files = [os.path.join(args.qtl_dir, file) for file in os.listdir(args.qtl_dir) if file.endswith('.txt')]
annotation_files = [os.path.join(args.annotation_dir, file) for file in os.listdir(args.annotation_dir) if file.endswith('.txt')]

# 读取所有文件并存储在字典中
qtl_data = {os.path.basename(file): pd.read_csv(file, header=None, names=['SNP']).squeeze("columns") for file in qtl_files}
annotation_data = {os.path.basename(file): pd.read_csv(file, header=None, names=['SNP']).squeeze("columns") for file in annotation_files}

# 结果存储列表
results = []

# 总SNP集合
all_snps = set(pd.read_csv(args.site_file, header=None).squeeze("columns"))

# 逐一计算每种QTL和功能注释的OR值及95%置信区间
for qtl_name, qtl_snps in qtl_data.items():
    qtl_snps_set = set(qtl_snps)
    num_var = len(qtl_snps_set)
    for annotation_name, annotation_snps in annotation_data.items():
        annotation_snps_set = set(annotation_snps)
        
        # 构建2x2列联表
        a = len(qtl_snps_set & annotation_snps_set)  # qtl组中与annot overlap的数量
        b = len(qtl_snps_set - annotation_snps_set)  # qtl组中与annot nonoverlap的数量
        c = len(annotation_snps_set - qtl_snps_set)  # nonqtl组中与annot overlap的数量
        d = len(all_snps - (qtl_snps_set | annotation_snps_set))  # nonqtl组中与annot nonoverlap的数量
        
        # 计算OR值和置信区间
        or_value, ci_low, ci_high, p_value = calculate_or_and_ci(a, b, c, d)
        
        # 存储结果
        results.append([qtl_name, annotation_name, or_value, ci_low, ci_high, p_value, num_var])

# 将结果转换为DataFrame并保存
results_df = pd.DataFrame(results, columns=['QTL', 'Annotation', 'OR', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_df.to_csv(f'{args.out_dir}/OR_results.csv', index=False)

print("计算完成，OR结果已保存到OR_results.csv文件中。")

# 结果存储列表
results = []

# 逐一计算每种QTL和功能注释的Fold enrichment及95%置信区间
for qtl_name, qtl_snps in qtl_data.items():
    qtl_snps_set = set(qtl_snps)
    num_var = len(qtl_snps_set)
    for annotation_name, annotation_snps in annotation_data.items():
        annotation_snps_set = set(annotation_snps)
        
        # 构建2x2列联表
        target_count = len(qtl_snps_set & annotation_snps_set)
        target_total = len(qtl_snps_set)
        background_count = len(annotation_snps_set)
        background_total = len(all_snps)
        
        # 计算FE值和置信区间
        fe_value, ci_low, ci_high, p_value = calculate_fold_enrichment(target_count, target_total, background_count, background_total)
        
        # 存储结果
        results.append([qtl_name, annotation_name, fe_value, ci_low, ci_high, p_value, num_var])

# 将结果转换为DataFrame并保存
results_df = pd.DataFrame(results, columns=['QTL', 'Annotation', 'Fold enrichment', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_df.to_csv(f'{args.out_dir}/FE_results.csv', index=False)

print("计算完成，Fold enrichment结果已保存到FE_results.csv文件中。")
