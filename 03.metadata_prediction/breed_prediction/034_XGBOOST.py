import os
import datatable as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import StandardScaler, LabelEncoder, label_binarize
import logging
import argparse
import xgboost as xgb

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s", datefmt="%I:%M:%S %p")

def main(encoded_csv_path, pheno_path, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Load encoded genotype data
    logging.info("Loading encoded genotype data with datatable...")
    dt_encoded = dt.fread(encoded_csv_path)
    df_encoded = dt_encoded.to_pandas()
    
    # Load phenotype data
    logging.info("Loading phenotype data with pandas...")
    pheno_data = pd.read_csv(pheno_path, sep='\t', index_col='SampleID')
    
    # Transpose and join genotype and phenotype data
    logging.info("Joining genotype and phenotype data...")
    df_encoded_transposed = df_encoded.set_index('Position').T
    df_joined = df_encoded_transposed.join(pheno_data, how='inner')
    
    if df_joined.empty:
        raise ValueError("The joined DataFrame is empty. Check that the 'SampleID' values in both the genotype and phenotype data match.")
    
    logging.info(f"Joined data shape: {df_joined.shape}")

    # Filter out phenotypes with less than 5 samples
    phenotype_counts = df_joined['Phenotype'].value_counts()
    valid_phenotypes = phenotype_counts[phenotype_counts > 5].index
    df_joined = df_joined[df_joined['Phenotype'].isin(valid_phenotypes)]
    
    logging.info(f"Filtered data shape after removing rare phenotypes: {df_joined.shape}")

    # Encode phenotype labels
    logging.info("Encoding phenotype labels...")
    label_encoder = LabelEncoder()
    df_joined['Phenotype'] = label_encoder.fit_transform(df_joined['Phenotype'])

    # Split dataset with stratification
    logging.info("Splitting dataset with stratification...")
    X = df_joined.drop('Phenotype', axis=1)
    y = df_joined['Phenotype']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    # Scale data
    logging.info("Scaling data...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Initialize XGBoost classifier
    logging.info("Initializing XGBoost classifier...")
    xgb_model = xgb.XGBClassifier(objective='binary:logistic', seed=42)

    # Hyperparameter tuning using GridSearchCV
    param_grid = {
        'n_estimators': [100, 200, 300],
        'learning_rate': [0.01, 0.1, 0.2],
        'max_depth': [3, 5, 7],
        'min_child_weight': [1, 3, 5],
        'subsample': [0.8, 0.9, 1.0],
        'colsample_bytree': [0.8, 0.9, 1.0]
    }

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    grid_search = GridSearchCV(estimator=xgb_model, param_grid=param_grid, scoring='accuracy', n_jobs=-1, cv=skf, verbose=2)
    grid_search.fit(X_train_scaled, y_train)

    logging.info("Best parameters found: %s", grid_search.best_params_)
    xgb_best = grid_search.best_estimator_

    # Fit the best model on the training data
    xgb_best.fit(X_train_scaled, y_train)

    # Predict and calculate accuracy on test data
    predictions = xgb_best.predict(X_test_scaled)
    accuracy = accuracy_score(y_test, predictions)
    logging.info(f"Model accuracy on test data: {accuracy:.4f}")

    logging.info("Generating classification report and confusion matrix...")
    class_report = classification_report(y_test, predictions, target_names=label_encoder.classes_)
    conf_matrix = confusion_matrix(y_test, predictions)

    logging.info(f"Classification Report:\n{class_report}")
    logging.info(f"Confusion Matrix:\n{conf_matrix}")

    evaluation_results_path = os.path.join(output_dir, 'evaluation_results.txt')
    with open(evaluation_results_path, 'w') as f:
        f.write(f"Model accuracy: {accuracy:.4f}\n\n")
        f.write(f"Classification Report:\n{class_report}\n")
        f.write(f"Confusion Matrix:\n{conf_matrix}\n")

    # Generate ROC and PR curves
    logging.info("Generating ROC and PR curves...")

    y_test_bin = label_binarize(y_test, classes=np.unique(y))
    y_pred_prob = xgb_best.predict_proba(X_test_scaled)
    n_classes = y_test_bin.shape[1]

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_pred_prob[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    fpr["micro"], tpr["micro"], _ = roc_curve(y_test_bin.ravel(), y_pred_prob.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    plt.figure()
    plt.plot(fpr["micro"], tpr["micro"], label='micro-average ROC curve (area = {0:0.2f})'.format(roc_auc["micro"]))
    for i in range(n_classes):
        plt.plot(fpr[i], tpr[i], lw=2, label='ROC curve of class {0} (area = {1:0.2f})'.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    roc_curve_path = os.path.join(output_dir, "roc_curve.png")
    plt.savefig(roc_curve_path)

    precision = dict()
    recall = dict()
    pr_auc = dict()
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(y_test_bin[:, i], y_pred_prob[:, i])
        pr_auc[i] = auc(recall[i], precision[i])

    precision["micro"], recall["micro"], _ = precision_recall_curve(y_test_bin.ravel(), y_pred_prob.ravel())
    pr_auc["micro"] = auc(recall["micro"], precision["micro"])

    plt.figure()
    plt.plot(recall["micro"], precision["micro"], label='micro-average PR curve (area = {0:0.2f})'.format(pr_auc["micro"]))
    for i in range(n_classes):
        plt.plot(recall[i], precision[i], lw=2, label='PR curve of class {0} (area = {1:0.2f})'.format(i, pr_auc[i]))

    plt.plot([0, 1], [1, 0], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall curve')
    plt.legend(loc="lower left")
    pr_curve_path = os.path.join(output_dir, "pr_curve.png")
    plt.savefig(pr_curve_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process encoded genotype and phenotype data for machine learning with XGBoost and 5-fold cross-validation.")
    parser.add_argument('encoded_csv_path', type=str, help='Path to the CSV file containing encoded genotypes.')
    parser.add_argument('pheno_path', type=str, help='Path to the phenotype CSV file.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output files.')
    args = parser.parse_args()

    main(args.encoded_csv_path, args.pheno_path, args.output_dir)

