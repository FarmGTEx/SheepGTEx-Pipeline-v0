import pandas as pd
import datatable as dt
import pickle
import argparse
import os
from sklearn.preprocessing import StandardScaler

def load_model(model_path):
    # 从文件加载训练好的模型
    with open(model_path, 'rb') as file:
        model = pickle.load(file)
    return model

def load_feature_names(feature_names_path):
    # 从文件加载特征名称
    with open(feature_names_path, 'rb') as file:
        feature_names = pickle.load(file)
    return feature_names

def load_scaler(scaler_path):
    # 从文件加载标准化器
    with open(scaler_path, 'rb') as file:
        scaler = pickle.load(file)
    return scaler

def process_data(encoded_csv_path, feature_names, scaler):
    # 读取编码后的数据并转换格式
    dt_encoded = dt.fread(encoded_csv_path)
    df_encoded = dt_encoded.to_pandas()
    # 确保数据的索引是位点信息
    df_encoded.set_index('Position', inplace=True)
    # 转置DataFrame以使样本成为行
    df_encoded_transposed = df_encoded.T
    # 确保特征顺序与训练时一致
    df_encoded_transposed = df_encoded_transposed.reindex(columns=feature_names)
    # 标准化特征，但先保存原始数据的索引
    original_index = df_encoded_transposed.index
    scaled_features = scaler.transform(df_encoded_transposed)
    return scaled_features, original_index  # 返回标准化特征和原始索引

def predict(model, data):
    # 使用模型对新数据进行预测
    return model.predict(data)

def main(model_path, feature_names_path, scaler_path, new_data_path, output_path):
    # 主函数
    model = load_model(model_path)
    feature_names = load_feature_names(feature_names_path)
    scaler = load_scaler(scaler_path)
    new_data, original_index = process_data(new_data_path, feature_names, scaler)  # 接收标准化数据和原始索引
    predictions = predict(model, new_data)
    # 保存预测结果
    predictions_path = os.path.join(output_path, 'predictions.csv')
    pd.Series(predictions, index=original_index).to_csv(predictions_path)  # 使用原始索引

if __name__ == '__main__':
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Use a trained SVM model to predict new encoded data.")
    parser.add_argument('model_path', type=str, help='Path to the trained model file.')
    parser.add_argument('feature_names_path', type=str, help='Path to the feature names file.')
    parser.add_argument('scaler_path', type=str, help='Path to the scaler file.')
    parser.add_argument('new_data_path', type=str, help='Path to the new encoded data CSV file.')
    parser.add_argument('output_path', type=str, help='Directory to save the prediction results.')
    args = parser.parse_args()

    main(args.model_path, args.feature_names_path, args.scaler_path, args.new_data_path, args.output_path)

