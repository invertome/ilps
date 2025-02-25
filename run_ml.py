#!/usr/bin/env python3
# Purpose: Train an ensemble ML model (RF + XGBoost) and predict ILP probabilities/novelty
# Inputs: 
#   - analysis/features.csv (candidate features)
#   - preprocess/ref_features.csv (reference features)
#   - preprocess/ref_labels.csv (reference labels)
#   - analysis/predictions.csv (output probabilities)
#   - analysis/novel_candidates.csv (output novelty flags)
#   - max CPUs (for parallel processing)
# Outputs: 
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
#   - output/shap_summary.png (SHAP plot for feature importance)
#   - output/rf_model.joblib, output/xgb_model.joblib (saved models)
# Notes: Uses grid search for hyperparameter tuning if models don’t exist; includes feature selection via SHAP; dynamic chunk size based on memory
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV
import shap
import matplotlib.pyplot as plt
import joblib
import sys
import os
import psutil

# Parse command-line arguments
candidate_features_file, ref_features_file, ref_labels_file, pred_file, novel_file, max_cpus = sys.argv[1:7]

# Dynamic chunk size based on available memory (approx. 10MB per chunk)
chunksize = int(psutil.virtual_memory().available / (1024 * 1024 * 10))

# Load datasets with chunked processing for memory efficiency
candidate_features = pd.concat([chunk for chunk in pd.read_csv(candidate_features_file, chunksize=chunksize)])
ref_features = pd.concat([chunk for chunk in pd.read_csv(ref_features_file, chunksize=chunksize)])
ref_labels = pd.read_csv(ref_labels_file, names=["label"])

# Prepare training and test datasets, dropping non-numeric columns
X_train_full = ref_features.drop(columns=["id", "domains", "motifs"])
y_train = ref_labels["label"]
X_test_full = candidate_features.drop(columns=["id", "domains", "motifs"])

# Define model file paths
rf_model_file = "output/rf_model.joblib"
xgb_model_file = "output/xgb_model.joblib"

# Load existing models or train new ones
if os.path.exists(rf_model_file) and os.path.exists(xgb_model_file):
    rf = joblib.load(rf_model_file)
    xgb = joblib.load(xgb_model_file)
else:
    # Random Forest hyperparameter tuning with full feature set
    rf_param_grid = {
        "n_estimators": [50, 100, 200],      # Number of trees
        "max_depth": [None, 10, 20],         # Maximum depth of trees
        "min_samples_split": [2, 5]          # Minimum samples to split a node
    }
    rf = GridSearchCV(RandomForestClassifier(n_jobs=int(max_cpus), random_state=42), rf_param_grid, cv=5, n_jobs=int(max_cpus))
    rf.fit(X_train_full, y_train)
    rf_full = rf.best_estimator_  # Keep full model for SHAP analysis

    # Feature selection using SHAP values
    explainer = shap.TreeExplainer(rf_full)
    shap_values = explainer.shap_values(X_train_full)
    shap_sum = pd.DataFrame(shap_values[1], columns=X_train_full.columns).abs().mean()
    top_features = shap_sum.nlargest(10).index  # Select top 10 features
    X_train = X_train_full[top_features]
    X_test = X_test_full[top_features]

    # Train optimized RF with selected features
    rf = GridSearchCV(RandomForestClassifier(n_jobs=int(max_cpus), random_state=42), rf_param_grid, cv=5, n_jobs=int(max_cpus))
    rf.fit(X_train, y_train)
    rf = rf.best_estimator_

    # XGBoost hyperparameter tuning with selected features
    xgb_param_grid = {
        "n_estimators": [50, 100, 200],      # Number of boosting rounds
        "max_depth": [3, 6, 9],              # Maximum tree depth
        "learning_rate": [0.01, 0.1]         # Step size shrinkage
    }
    xgb = GridSearchCV(XGBClassifier(n_jobs=int(max_cpus), random_state=42), xgb_param_grid, cv=5, n_jobs=int(max_cpus))
    xgb.fit(X_train, y_train)
    xgb = xgb.best_estimator_

    # Save trained models for reuse
    os.makedirs("output", exist_ok=True)
    joblib.dump(rf, rf_model_file)
    joblib.dump(xgb, xgb_model_file)

    # Generate SHAP summary plot with full model for visualization
    shap.summary_plot(shap_values, X_train_full, show=False)
    plt.savefig("output/shap_summary.png", bbox_inches="tight")
    plt.close()

# Predict probabilities with ensemble (average of RF and XGBoost)
X_test = X_test_full[top_features]  # Use selected features for prediction
preds = (rf.predict_proba(X_test)[:, 1] + xgb.predict_proba(X_test)[:, 1]) / 2

# Identify novel ILPs (high probability, low similarity to known ILPs)
novel = ((preds > 0.7) & (candidate_features["hhsearch_prob"] < 70) & (candidate_features["blast_identity"] < 30)).astype(int)

# Write predictions and novelty flags to CSV
pd.Series(preds).to_csv(pred_file, index=False)
pd.Series(novel).to_csv(novel_file, index=False)
