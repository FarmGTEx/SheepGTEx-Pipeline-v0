import pandas as pd
from scipy.stats import fisher_exact
import numpy as np
import os
#import statsmodels.api as sm
import argparse


# 定义计算OR值和95%置信区间的函数
def calculate_or_and_ci(a, b, c, d):
    # 检查零值
    if a == 0 or b == 0 or c == 0 or d == 0:
        or_value = float('inf') if (a * d > b * c) else 0
        ci_low, ci_high = (0, float('inf'))
        p_value = 1.0 # 如果不能计算，则P值设为1

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
        ci_low = np.exp(np.log(or_value) - 1.96 * se)
        ci_high = np.exp(np.log(or_value) + 1.96 * se)
    
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
        se = np.sqrt((1 - target_ratio) / target_count + (1 - background_ratio) / background_count)
        ci_low = np.exp(np.log(fe_value) - 1.96 * se)
        ci_high = np.exp(np.log(fe_value) + 1.96 * se)
    
    return fe_value, ci_low, ci_high, p_value

# 获得control SNP
def get_control_snps(maf_ld_data, ld_sd, qtl_snps_set, maf_match, ld_match, permutations):
    qtl_snps_set_control = set()
    for snp in list(qtl_snps_set):
        # 找到该 SNP 的信息
        if snp not in maf_ld_data.index:
            continue
        focal_maf = maf_ld_data.at[snp, "MAF"]
        focal_ld = maf_ld_data.at[snp, "ldscore"]
        maf_low = focal_maf - maf_match
        maf_high = focal_maf + maf_match
        ld_low = focal_ld - ld_match * ld_sd
        ld_high = focal_ld + ld_match * ld_sd
        
        # 筛选符合条件的 SNP
        candidates = maf_ld_data[
            (maf_ld_data["MAF"].between(maf_low, maf_high)) &
            (maf_ld_data["ldscore"].between(ld_low, ld_high))
        ]
        candidates = candidates.drop(index=snp, errors="ignore")
        candidates = candidates[~candidates.index.isin(qtl_snps_set)]

        # 随机选择n个
        if len(candidates) >= permutations:
            matched = np.random.choice(candidates.index, permutations, replace=False)
        else:
            matched = candidates.index.values
        
        qtl_snps_set_control.update(matched)
    return qtl_snps_set_control

    
# 文件目录设置
parser = argparse.ArgumentParser(description='calculate the OR/FE of molQTL.')
parser.add_argument('qtl_dir', help='directory of files (end with .txt) of molQTL. variant_id should be included wihout headers')
parser.add_argument('annotation_dir', help='directory of files (end with .txt) of functional annotation sites.')
parser.add_argument('site_file', help='file of all sites. variant_id should be included wihout headers')
parser.add_argument('out_dir', help='directory to save results (end with .csv).')
parser.add_argument('--maf_ld_file', help='file of ldscore generated from GCTA (end with .score.ld). At lead three fields should be included with header: SNP, MAF and ldscore')
parser.add_argument('--maf_match', type=float, default=0.02, help='control variants within the range of MAF of the target variants. (default: 0.02)')
parser.add_argument('--ld_match', type=float, default=0.1, help='control variants with LD score within the s.d. of the focal variant’s LD score. (default: 0.1)')
args = parser.parse_args()

# 获取目录中的文件列表
qtl_files = [os.path.join(args.qtl_dir, file) for file in os.listdir(args.qtl_dir) if file.endswith('.txt')]
annotation_files = [os.path.join(args.annotation_dir, file) for file in os.listdir(args.annotation_dir) if file.endswith('.txt')]

# 读取所有文件并存储在字典中
qtl_data = {os.path.basename(file): pd.read_csv(file, header=None, names=['SNP']).squeeze("columns") for file in qtl_files}
annotation_data = {os.path.basename(file): pd.read_csv(file, header=None, names=['SNP']).squeeze("columns") for file in annotation_files}
if args.maf_ld_file:
    maf_ld_data = pd.read_csv(args.maf_ld_file, sep=' ').dropna().set_index("SNP")
    ld_sd = maf_ld_data['ldscore'].std()

# 总SNP集合
all_snps = set(pd.read_csv(args.site_file, header=None).squeeze("columns"))

# 结果存储列表
results_or = []
results_fe = []
results_or_control = []
results_fe_control = []
qtl_snps_set_control = set()

# 逐一计算每种QTL和功能注释的OR/FE值及95%置信区间
for qtl_name, qtl_snps in qtl_data.items():
    qtl_snps_set = set(qtl_snps)
    
    # 提取MAF-matched和LD-matched control SNPs，保证至少10个重复，且大于1000个位点
    if args.maf_ld_file:
        num_var = len(qtl_snps_set)
        if num_var > 10000:
            # downsampling to 10000 for focal SNPs
            qtl_snps_set_control = get_control_snps(maf_ld_data, ld_sd, set(np.random.choice(list(qtl_snps_set), 10000, replace=False)), args.maf_match, args.ld_match, 10)
        elif num_var > 1000:
            qtl_snps_set_control = get_control_snps(maf_ld_data, ld_sd, qtl_snps_set, args.maf_match, args.ld_match, 10)
        elif num_var > 10:
            qtl_snps_set_control = get_control_snps(maf_ld_data, ld_sd, qtl_snps_set, args.maf_match, args.ld_match, 100)
        else:
            qtl_snps_set_control = get_control_snps(maf_ld_data, ld_sd, qtl_snps_set, args.maf_match, args.ld_match, 1000)
        print(qtl_name, num_var, len(qtl_snps_set_control))
    
    for annotation_name, annotation_snps in annotation_data.items():
        annotation_snps_set = set(annotation_snps)
        num_var = len(qtl_snps_set)
        
        # OR
        ## 构建2x2列联表
        a = len(qtl_snps_set & annotation_snps_set)  # qtl组中与annot overlap的数量
        b = len(qtl_snps_set - annotation_snps_set)  # qtl组中与annot nonoverlap的数量
        c = len(annotation_snps_set - qtl_snps_set)  # nonqtl组中与annot overlap的数量
        d = len(all_snps - (qtl_snps_set | annotation_snps_set))  # nonqtl组中与annot nonoverlap的数量
        ## 计算OR值和置信区间
        or_value, ci_low, ci_high, p_value = calculate_or_and_ci(a, b, c, d)
        ## 存储结果
        results_or.append([qtl_name, annotation_name, or_value, ci_low, ci_high, p_value, num_var])
        
        # FE
        ## 构建2x2列联表
        target_count = len(qtl_snps_set & annotation_snps_set)
        target_total = len(qtl_snps_set)
        background_count = len(annotation_snps_set)
        background_total = len(all_snps)
        ## 计算FE值和置信区间
        fe_value, ci_low, ci_high, p_value = calculate_fold_enrichment(target_count, target_total, background_count, background_total)
        ## 存储结果
        results_fe.append([qtl_name, annotation_name, fe_value, ci_low, ci_high, p_value, num_var])

        if qtl_snps_set_control:
            num_var = len(qtl_snps_set_control)
            
            # OR (control)
            ## 构建2x2列联表
            a = len(qtl_snps_set_control & annotation_snps_set)  # qtl组中与annot overlap的数量
            b = len(qtl_snps_set_control - annotation_snps_set)  # qtl组中与annot nonoverlap的数量
            c = len(annotation_snps_set - qtl_snps_set_control)  # nonqtl组中与annot overlap的数量
            d = len(all_snps - (qtl_snps_set_control | annotation_snps_set))  # nonqtl组中与annot nonoverlap的数量
            ## 计算OR值和置信区间
            or_value, ci_low, ci_high, p_value = calculate_or_and_ci(a, b, c, d)
            ## 存储结果
            results_or_control.append([qtl_name, annotation_name, or_value, ci_low, ci_high, p_value, num_var])
            
            # FE (control)
            ## 构建2x2列联表
            target_count = len(qtl_snps_set_control & annotation_snps_set)
            target_total = len(qtl_snps_set_control)
            background_count = len(annotation_snps_set)
            background_total = len(all_snps)
            ## 计算FE值和置信区间
            fe_value, ci_low, ci_high, p_value = calculate_fold_enrichment(target_count, target_total, background_count, background_total)
            ## 存储结果
            results_fe_control.append([qtl_name, annotation_name, fe_value, ci_low, ci_high, p_value, num_var])

# 将结果转换为DataFrame并保存
results_or_df = pd.DataFrame(results_or, columns=['QTL', 'Annotation', 'Odds ratio', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_or_df.to_csv(f'{args.out_dir}/OR_results.csv', index=False)

results_fe_df = pd.DataFrame(results_fe, columns=['QTL', 'Annotation', 'Fold enrichment', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_fe_df.to_csv(f'{args.out_dir}/FE_results.csv', index=False)

results_or_control_df = pd.DataFrame(results_or_control, columns=['QTL', 'Annotation', 'Odds ratio', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_or_control_df.to_csv(f'{args.out_dir}/OR_control_results.csv', index=False)

results_fe_control_df = pd.DataFrame(results_fe_control, columns=['QTL', 'Annotation', 'Fold enrichment', '95% CI Low', '95% CI High', 'P-value', 'num_var'])
results_fe_control_df.to_csv(f'{args.out_dir}/FE_control_results.csv', index=False)
