from joblib import load
import pandas as pd
import sys

# 确保至少提供了一个命令行参数（除了脚本名称本身）
if len(sys.argv) < 2:
    print("Usage: python script.py <filename_part>")
    sys.exit(1)  # 退出程序，并指示错误
#获取命令行参数
filename_part = sys.argv[1]

age_without = f"{filename_part}_without_age.txt"
age_predict=f"{filename_part}_without_predict_age.txt"
# 模型文件路径
model_path = f"{filename_part}_random_forest.joblib"
# 使用 joblib 加载模型
rf_model = load(model_path)

# 加载缺少age的数据
data_without_age = pd.read_csv(age_without, header=0,sep='\t')

# 假设你知道哪些列是特征列
# 这里需要保证这些列与训练模型时使用的特征列相匹配
X_new = data_without_age.drop(columns=['sample', 'age'])  # 如果sample列不是特征，需要删除

# 使用模型进行预测
y_new_predictions = rf_model.predict(X_new)

# 保存预测结果到文件
with open(age_predict, 'w') as f:
    for prediction in y_new_predictions:
        f.write(f"{prediction}\n")

