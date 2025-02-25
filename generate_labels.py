#!/usr/bin/env python3
# Purpose: Generate binary labels for reference sequences (1 for ILP, 0 for non-ILP)
# Inputs: 
#   - input/ref_ILPs.fasta (annotated sequences)
#   - preprocess/ref_features.csv (feature matrix)
# Output: preprocess/ref_labels.csv (labels)
# Notes: Matches sequence IDs to ensure alignment with features
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import pandas as pd
import sys

ref_fasta, features_file, output_labels = sys.argv[1:4]
features = pd.read_csv(features_file)
labels = []

# Assign labels based on [ILP] tag in description
for rec in SeqIO.parse(ref_fasta, "fasta"):
    label = 1 if "[ILP]" in rec.description else 0
    labels.append({"id": rec.id, "label": label})

# Merge with feature IDs and fill missing with 0 (non-ILP)
labels_df = pd.DataFrame(labels)
labels_df = labels_df.merge(features[["id"]], on="id", how="right").fillna(0)
labels_df["label"].to_csv(output_labels, index=False)
