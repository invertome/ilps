#!/bin/bash
# Purpose: Annotate candidate sequences with ML predictions and identify novel ILPs
# Inputs: 
#   - analysis/all_candidates.fasta (combined candidate sequences)
#   - preprocess/ref_features.csv (reference features)
#   - preprocess/ref_labels.csv (reference labels)
# Outputs: 
#   - analysis/features.csv (candidate features)
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
# Log: pipeline.log
# Notes: Uses pre-trained ML models for efficiency; prepares data for phylogenetics
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f analysis/.done_annotate ]; then
    echo "$(date) - Skipping 04_annotate_and_novel.sh (already done)" >> pipeline.log
    exit 0
fi

mkdir -p analysis
max_cpus=$(nproc)
echo "$(date) - Starting 04_annotate_and_novel.sh" >> pipeline.log

# Combine all candidate sequences into a single file
cat candidates/[0-9]*_candidates.fasta > analysis/all_candidates.fasta || { echo "$(date) - cat failed" >> pipeline.log; exit 1; }

# Extract features for ML prediction
python extract_features.py analysis/all_candidates.fasta > analysis/features.csv || { echo "$(date) - extract_features.py failed" >> pipeline.log; exit 1; }

# Run ML model to predict ILP probabilities and novelty
python run_ml.py analysis/features.csv preprocess/ref_features.csv preprocess/ref_labels.csv analysis/predictions.csv analysis/novel_candidates.csv "$max_cpus" || { echo "$(date) - run_ml.py failed" >> pipeline.log; exit 1; }

echo "$(date) - 04_annotate_and_novel.sh completed" >> pipeline.log
touch analysis/.done_annotate
