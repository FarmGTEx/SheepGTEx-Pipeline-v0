
import datatable as dt
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
import logging
import argparse
import os
import pickle

# Setting up logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s", datefmt="%I:%M:%S %p")

def main(encoded_csv_path, pheno_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Loading encoded genotype data with datatable...")
    dt_encoded = dt.fread(encoded_csv_path)
    df_encoded = dt_encoded.to_pandas()

    logging.info("Loading phenotype data with pandas...")
    pheno_data = pd.read_csv(pheno_path, sep='\t', index_col='SampleID')

    logging.info("Joining genotype and phenotype data...")
    df_encoded_transposed = df_encoded.set_index('Position').T
    df_joined = df_encoded_transposed.join(pheno_data, how='inner')

    if df_joined.empty:
        raise ValueError("The joined DataFrame is empty. Check that the 'SampleID' values in both the genotype and phenotype data match.")

    logging.info(f"Joined data shape: {df_joined.shape}")

    logging.info("Splitting dataset...")
    X = df_joined.drop('Phenotype', axis=1)
    y = df_joined['Phenotype']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    logging.info("Scaling data...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    logging.info("Initializing KNN classifier and grid search...")
    knn = KNeighborsClassifier()
    param_grid = {
        'n_neighbors': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21],
        'weights': ['uniform', 'distance'],
        'metric': ['euclidean', 'manhattan']
    }
    grid_search = GridSearchCV(knn, param_grid, cv=5, n_jobs=-1, verbose=2)
    grid_search.fit(X_train_scaled, y_train)

    logging.info("Saving the best estimator...")
    best_estimator_path = os.path.join(output_dir, 'best_estimator.pkl')
    pd.to_pickle(grid_search.best_estimator_, best_estimator_path)

    logging.info("Saving feature names...")
    feature_names_path = os.path.join(output_dir, 'feature_names.pkl')
    with open(feature_names_path, 'wb') as f:
        pickle.dump(X_train.columns, f)

    logging.info("Predicting with the best estimator...")
    predictions = grid_search.best_estimator_.predict(X_test_scaled)
    predictions_path = os.path.join(output_dir, 'predictions.csv')
    pd.Series(predictions, index=X_test.index).to_csv(predictions_path)

    accuracy = accuracy_score(y_test, predictions)
    logging.info(f"Model accuracy: {accuracy:.4f}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process encoded genotype and phenotype data for machine learning.")
    parser.add_argument('encoded_csv_path', type=str, help='Path to the CSV file containing encoded genotypes.')
    parser.add_argument('pheno_path', type=str, help='Path to the phenotype CSV file.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output files.')
    args = parser.parse_args()

    main(args.encoded_csv_path, args.pheno_path, args.output_dir)

