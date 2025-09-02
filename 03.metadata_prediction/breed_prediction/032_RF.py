import datatable as dt
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import logging
import argparse
import os
import pickle

# 设置日志记录的配置
logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s", datefmt="%I:%M:%S %p")

def main(encoded_csv_path, pheno_path, output_dir):
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 使用datatable读取编码后的基因型数据
    logging.info("Loading encoded genotype data with datatable...")
    dt_encoded = dt.fread(encoded_csv_path)
    df_encoded = dt_encoded.to_pandas()

    # 使用pandas读取表型数据
    logging.info("Loading phenotype data with pandas...")
    pheno_data = pd.read_csv(pheno_path, sep='\t', index_col='SampleID')

    # 基因型数据转置，使样本ID成为行索引
    df_encoded_transposed = df_encoded.set_index('Position').T

    # 结合转置后的基因型数据和表型数据
    logging.info("Joining genotype and phenotype data...")
    df_joined = df_encoded_transposed.join(pheno_data, how='inner')

    if df_joined.empty:
        raise ValueError("The joined DataFrame is empty. Check that the 'SampleID' values in both the genotype and phenotype data match.")

    logging.info(f"Joined data shape: {df_joined.shape}")

    # 准备数据集
    logging.info("Preparing dataset for cross-validation...")
    X = df_joined.drop('Phenotype', axis=1)
    y = df_joined['Phenotype']

    # 数据标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # 初始化随机森林分类器
    rf = RandomForestClassifier(random_state=42)

    # 进行5折交叉验证
    logging.info("Performing 5-fold cross-validation...")
    cv_scores = cross_val_score(rf, X_scaled, y, cv=5)
    logging.info(f"Cross-validation scores: {cv_scores}")

    # 计算平均准确度
    mean_accuracy = cv_scores.mean()
    logging.info(f"Mean accuracy from cross-validation: {mean_accuracy:.4f}")

    # 训练模型并保存
    rf.fit(X_scaled, y)
    model_path = os.path.join(output_dir, 'random_forest_model.pkl')
    with open(model_path, 'wb') as f:
        pickle.dump(rf, f)
    logging.info(f"Trained model saved to {model_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process encoded genotype and phenotype data for machine learning with Random Forest.")
    parser.add_argument('encoded_csv_path', type=str, help='Path to the CSV file containing encoded genotypes.')
    parser.add_argument('pheno_path', type=str, help='Path to the phenotype CSV file.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output files.')
    args = parser.parse_args()

    main(args.encoded_csv_path, args.pheno_path, args.output_dir)

