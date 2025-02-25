#!/usr/bin/env python3
# Purpose: Calculate minimum and maximum lengths of reference ILP prepropeptides for filtering
# Inputs: input/ref_ILPs.fasta (annotated ILP/non-ILP sequences)
# Output: preprocess/ref_lengths.txt (min and max lengths as two lines)
# Notes: Used to set dynamic length thresholds in preprocessing step
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import sys

ref_fasta, output_file = sys.argv[1:3]
# Filter for ILP sequences and calculate lengths
lengths = [len(rec.seq) for rec in SeqIO.parse(ref_fasta, "fasta") if "[ILP]" in rec.description]
# Write min and max lengths to file
with open(output_file, "w") as f:
    f.write(f"{min(lengths)}\n{max(lengths)}")
