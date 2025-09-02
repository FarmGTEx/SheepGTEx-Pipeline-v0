# -*- coding:utf-8 -*-
# pandas 完成数据清洗和预备
# 使用随机森林方法进行数据训练和预测以及结果精度分析
import sys

# 确保至少提供了一个命令行参数（除了脚本名称本身）
if len(sys.argv) < 2:
    print("Usage: python script.py <filename_part>")
    sys.exit(1)  # 退出程序，并指示错误
# 获取命令行参数
filename_part = sys.argv[1]
inputfile=f"{filename_part}.with_age.txt"
outputfile=f"{filename_part}.out.txt"
outputpng=f"{filename_part}.png"
outputpredict=f"{filename_part}.predict.txt"
age_with=f"{filename_part}.middle.age_with.txt"
age_without = f"{filename_part}_without_age.txt"
random=f"{filename_part}_random_forest.joblib"
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn import metrics
from sklearn.model_selection import train_test_split
# 该函数用与将已有的数据进行分离，区分开训练集和测试集
from sklearn.preprocessing import label_binarize
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier

if __name__ == '__main__':
    pd.set_option('display.width', 300)
    pd.set_option('display.max_columns', 300)
    # 设置numpy输出的参数宽度和最大的展示的个数 可以不用管

    data = pd.read_csv(inputfile, header=0)

    # 读取已有的数据
    n_columns = len(data.columns)

    # 筛选出含有缺失值的行
    data_age_na = data[data['age'].isna()]
    data_age_na.to_csv(age_without, sep='\t', index=False)
    ##删除缺失age的行
    data.dropna(subset=['age'], inplace=True)
    print(data)

data.to_csv(age_with, sep='\t', index=False)

# 首先，确保目标变量和特征集被正确地分开
# 假设最后一列是你需要预测的内容
X = data.drop(columns=['sample', 'age'])  # 删除非特征列
y = data['age']  # 目标变量

# 划分数据集为训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 创建随机森林模型
rf_model = RandomForestClassifier(n_estimators=1500,      # 森林中的树木数量
                                  min_samples_leaf=3,     # 叶节点上的最小样本数
                                  max_features='sqrt',    # 查找最佳分割时要考虑的特征数量
                                  max_depth=48,           # 树的最大深度
                                  bootstrap=False,         # 是否使用bootstrap样本
                                 n_jobs=-1)

# 训练模型
rf_model.fit(X_train, y_train)

# 进行预测
y_pred = rf_model.predict(X_test)

# 评估模型精度
accuracy = metrics.accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy}")

# 如果你感兴趣，还可以打印更多的性能指标
print(metrics.classification_report(y_test, y_pred))



from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# 保存结果
accuracy = accuracy_score(y_test, y_pred)
report = classification_report(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

with open(outputfile, 'w') as f:
    f.write(f"Accuracy: {accuracy}\n")
    f.write("\nClassification Report:\n")
    f.write(report)
    f.write("\nConfusion Matrix:\n")
    f.write(str(conf_matrix))


from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

# 假设 y_test 和 y_pred 已经被定义
cm = confusion_matrix(y_test, y_pred)

plt.figure(figsize=(10, 8))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel('Predicted')
plt.ylabel('True')

# 保存热力图为PNG图像文件
plt.savefig(outputpng)

# 将预测结果保存到文件中
with open(outputpredict, 'w') as f:
    for pred in y_pred:
        f.write(str(pred) + '\n')

from joblib import dump, load

# 训练模型后保存
dump(rf_model, random) 

