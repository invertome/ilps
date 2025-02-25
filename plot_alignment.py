#!/usr/bin/env python3
# Purpose: Generate sequence logo plots for clade alignments
# Inputs: 
#   - Aligned FASTA (e.g., analysis/aligned_prepro/ete_clade_0_aligned.fasta)
#   - Output plot path (e.g., analysis/plots_prepro/ete_clade_0_logo.png)
# Output: Sequence logo PNG file
# Notes: Uses logomaker for visualization
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import sys

input_fasta, output_plot = sys.argv[1:3]
seqs = [str(rec.seq) for rec in SeqIO.parse(input_fasta, "fasta")]
if not seqs:
    print(f"No sequences in {input_fasta}")
    sys.exit(0)

# Create frequency matrix for logo
matrix = pd.DataFrame([list(seq) for seq in seqs]).apply(pd.value_counts).fillna(0).T
logo = logomaker.Logo(matrix, color_scheme="chemistry")
plt.title(f"Sequence Logo for {os.path.basename(input_fasta)}")
plt.savefig(output_plot, bbox_inches="tight")
plt.close()
