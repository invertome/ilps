#!/usr/bin/env python3
# Purpose: Train ensemble ML model (RF + XGBoost) with cross-validation and predict ILP probabilities/novelty
# Inputs: analysis/features.csv, preprocess/ref_features.csv, preprocess/ref_labels.csv, output paths, max CPUs
# Outputs: analysis/predictions.csv, analysis/novel_candidates.csv, output/shap_summary.png, output/rf_model.joblib, output/xgb_model.joblib, output/ml_metrics.txt
# Notes: Includes 5-fold CV and metrics (AUC, precision, recall)
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.metrics import roc_auc_score, precision_score, recall_score
import shap
import matplotlib.pyplot as plt
import joblib
import sys
import os
import psutil
import numpy as np
import yaml

# Parse arguments and config
candidate_features_file, ref_features_file, ref_labels_file, pred_file, novel_file, max_cpus = sys.argv[1:7]
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Dynamic chunk size
chunksize = int(psutil.virtual_memory().available / (1024 * 1024 * 10))

# Load datasets
candidate_features = pd.concat([chunk for chunk in pd.read_csv(candidate_features_file, chunksize=chunksize)])
ref_features = pd.concat([chunk for chunk in pd.read_csv(ref_features_file, chunksize=chunksize)])
ref_labels = pd.read_csv(ref_labels_file, names=["label"])

X_train_full = ref_features.drop(columns=["id", "domains", "motifs"])
y_train = ref_labels["label"]
X_test_full = candidate_features.drop(columns=["id", "domains", "motifs"])

# Model file paths
rf_model_file = "output/rf_model.joblib"
xgb_model_file = "output/xgb_model.joblib"
metrics_file = "output/ml_metrics.txt"

if os.path.exists(rf_model_file) and os.path.exists(xgb_model_file):
    rf = joblib.load(rf_model_file)
    xgb = joblib.load(xgb_model_file)
else:
    # Hyperparameter grids
    rf_param_grid = {"n_estimators": [50, 100, 200], "max_depth": [None, 10, 20], "min_samples_split": [2, 5]}
    xgb_param_grid = {"n_estimators": [50, 100, 200], "max_depth": [3, 6, 9], "learning_rate": [0.01, 0.1]}

    # Initial RF training for SHAP
    rf = GridSearchCV(RandomForestClassifier(n_jobs=int(max_cpus), random_state=42), rf_param_grid, cv=5, n_jobs=int(max_cpus))
    rf.fit(X_train_full, y_train)
    rf_full = rf.best_estimator_

    # SHAP feature selection
    explainer = shap.TreeExplainer(rf_full)
    shap_values = explainer.shap_values(X_train_full)
    shap_sum = pd.DataFrame(shap_values[1], columns=X_train_full.columns).abs().mean()
    top_features = shap_sum.nlargest(10).index
    X_train = X_train_full[top_features]
    X_test = X_test_full[top_features]

    # Optimized RF with CV
    rf = GridSearchCV(RandomForestClassifier(n_jobs=int(max_cpus), random_state=42), rf_param_grid, cv=5, n_jobs=int(max_cpus))
    rf.fit(X_train, y_train)
    rf = rf.best_estimator_
    rf_cv_auc = cross_val_score(rf, X_train, y_train, cv=5, scoring="roc_auc").mean()
    rf_cv_precision = cross_val_score(rf, X_train, y_train, cv=5, scoring="precision").mean()
    rf_cv_recall = cross_val_score(rf, X_train, y_train, cv=5, scoring="recall").mean()

    # XGBoost with CV
    xgb = GridSearchCV(XGBClassifier(n_jobs=int(max_cpus), random_state=42), xgb_param_grid, cv=5, n_jobs=int(max_cpus))
    xgb.fit(X_train, y_train)
    xgb = xgb.best_estimator_
    xgb_cv_auc = cross_val_score(xgb, X_train, y_train, cv=5, scoring="roc_auc").mean()
    xgb_cv_precision = cross_val_score(xgb, X_train, y_train, cv=5, scoring="precision").mean()
    xgb_cv_recall = cross_val_score(xgb, X_train, y_train, cv=5, scoring="recall").mean()

    # Save models
    os.makedirs("output", exist_ok=True)
    joblib.dump(rf, rf_model_file)
    joblib.dump(xgb, xgb_model_file)

    # SHAP plot
    shap.summary_plot(shap_values, X_train_full, show=False)
    plt.savefig("output/shap_summary.png", bbox_inches="tight")
    plt.close()

    # Write metrics
    with open(metrics_file, "w") as f:
        f.write(f"RF CV AUC: {rf_cv_auc:.3f}\nRF CV Precision: {rf_cv_precision:.3f}\nRF CV Recall: {rf_cv_recall:.3f}\n")
        f.write(f"XGB CV AUC: {xgb_cv_auc:.3f}\nXGB CV Precision: {xgb_cv_precision:.3f}\nXGB CV Recall: {xgb_cv_recall:.3f}\n")

# Predict probabilities
X_test = X_test_full[top_features]
preds = (rf.predict_proba(X_test)[:, 1] + xgb.predict_proba(X_test)[:, 1]) / 2

# Identify novel ILPs
novel = ((preds > config["ilp_prob_threshold"]) & 
         (candidate_features["hhsearch_prob"] < config["hhsearch_novel_threshold"]) & 
         (candidate_features["blast_identity"] < config["blast_novel_threshold"])).astype(int)

pd.Series(preds).to_csv(pred_file, index=False)
pd.Series(novel).to_csv(novel_file, index=False)
