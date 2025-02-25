#!/usr/bin/env python3
# Purpose: Filter candidate sequences to ILPs based on ML probability (> 0.7)
# Inputs: 
#   - analysis/all_candidates.fasta (combined candidates)
#   - analysis/predictions.csv (ML probabilities)
# Output: analysis/ilp_candidates.fasta (filtered ILP sequences)
# Notes: Ensures only high-confidence ILPs proceed to phylogenetics
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import pandas as pd
import sys

input_fasta, pred_file, output_fasta = sys.argv[1:4]
preds = pd.read_csv(pred_file, names=["ML_Prob"])
seqs = list(SeqIO.parse(input_fasta, "fasta"))
# Filter sequences with ML probability > 0.7
ilp_seqs = [seq for seq, pred in zip(seqs, preds["ML_Prob"]) if pred > 0.7]
SeqIO.write(ilp_seqs, output_fasta, "fasta")
