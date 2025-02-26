#!/bin/bash
# Purpose: Step 03 - Annotate candidates with ML predictions and identify novel ILPs (initial pass)
# Inputs: candidates/[0-9]*_candidates.fasta, preprocess/ref_features.csv, preprocess/ref_labels.csv
# Outputs: analysis/all_candidates.fasta, analysis/predictions.csv, analysis/novel_candidates.csv
# Config: config.yaml (max_cpus)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Initial ML pass without structural features; structural annotation follows in 04
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 03_annotate_and_novel.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f analysis/.done_annotate ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 03_annotate_and_novel.sh (already done)" >> pipeline.log
    exit 0
fi

# Check input files
for file in candidates/[0-9]*_candidates.fasta preprocess/ref_features.csv preprocess/ref_labels.csv; do
    if ! ls "$file" >/dev/null 2>&1; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: $file not found" >> pipeline.log
        exit 1
    fi
done

mkdir -p analysis

# Combine all candidate FASTA files
cat candidates/[0-9]*_candidates.fasta > analysis/all_candidates.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to combine candidates" >> pipeline.log; exit 1; }

# Extract initial sequence-based features (no PDBs yet)
python extract_features.py analysis/all_candidates.fasta analysis/features_initial.csv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - extract_features.py failed" >> pipeline.log; exit 1; }

# Run ML with sequence-based features only
python run_ml.py analysis/features_initial.csv preprocess/ref_features.csv preprocess/ref_labels.csv analysis/predictions.csv analysis/novel_candidates.csv "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - run_ml.py failed" >> pipeline.log; exit 1; }

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 03_annotate_and_novel.sh completed in ${runtime}s" >> pipeline.log
touch analysis/.done_annotate
