import os
import random
import pandas as pd
import gzip
import sys
import io
import numpy as np

#-------------------------------------
#--- parameters
# work path
# sample size

# group CV
group_method = "k_fold"
num_splits = 10
num_repeats = 5
# feature selection
feature_method = sys.argv[1]
num_features = int(sys.argv[2])
# classification
class_label = "Male"
class_model = sys.argv[3]

#-------------------------------------
#--- define the work path
wkdir = "./"
os.chdir(wkdir)
output_path = wkdir
sys.path.append(wkdir)
#--- define the path
import classify_function

input_file = sys.argv[4]

output_file = "{0}.F{1}.{2}.csv".format(feature_method, num_features, class_model)
featur_file = "{0}.F{1}.{2}.feature.csv".format(feature_method, num_features, class_model)
#--- read gene expression data
data = pd.read_csv(input_file, delimiter='\t', header=0, index_col=0)

# random.seed(20240321)
# data = data.sample(n=1000)

#--- divide the label (y) and features (X)
X = data.drop("Gender", axis = 1)
y = data["Gender"]
# log2(expr + 1) transfer 
X = np.log2(X + 1)
#--- transfer from data.frame to array
X = X.values
y = y.values

#--- save the prediction results
df_titles = pd.DataFrame({
    'Fold': [],
    'Accuracy': [],
    'Precision': [],
    'Recall': [],
    'F1': [],
    'feature_select_time': [],
    'predict_time': []
})
df_titles.to_csv("{0}/{1}".format(output_path, output_file), index = False)

#--- save the feature information
f_titles = pd.DataFrame({
        'Fold': [],
        'Feature_index': []
    })
f_titles.to_csv("{0}/{1}".format(output_path, featur_file), index = False)

# data grouped
data_splits = classify_function.split_data(X, y, method = group_method, num_splits = num_splits, num_repeats = num_repeats)

for i, (X_train, X_test, y_train, y_test) in enumerate(data_splits):
    print(f"Fold {i+1}:")
    #--- feature selection
    X_train_selected, X_test_selected, selected_feature_indices, feature_select_time = classify_function.select_features(X_train, y_train, X_test, method = feature_method, num_features = num_features)
    #--- clasification
    accuracy, cm, report, precision, recall, f1, predict_time = classify_function.predict(X_train_selected, y_train, X_test_selected, y_test, label = class_label, model = class_model)
    #--- save the results
    df_i = pd.DataFrame({
        'Fold': [i+1],
        'Accuracy': [accuracy],
        'Precision': [precision],
        'Recall': [recall],
        'F1': [f1],
        'feature_select_time': [feature_select_time],
        'predict_time': [predict_time]
    })
    df_i.to_csv("{0}/{1}".format(output_path, output_file), mode='a', header = False, index = False)

    f_i = pd.DataFrame({
        'Fold': [i + 1],
        'Feature_index': [selected_feature_indices]
    })
    f_i.to_csv("{0}/{1}".format(output_path, featur_file), mode='a', header = False, index = False)
