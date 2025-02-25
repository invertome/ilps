#!/usr/bin/env python3
# Purpose: Identify clades from phylogenetic trees for motif analysis
# Inputs: 
#   - Tree file (e.g., analysis/prepro_alignment.fasta.treefile)
#   - Clades directory (e.g., clades_ete_prepro/)
#   - Unaligned directory (e.g., analysis/unaligned_prepro/)
#   - Aligned directory (e.g., analysis/aligned_prepro/)
# Outputs: 
#   - clades_dir/clade_*.fasta (clade sequences)
#   - unaligned_dir/clade_*.fasta (unaligned clade sequences)
# Notes: Uses ETE3 with a support threshold of 70
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from ete3 import Tree
from Bio import SeqIO
import sys
import os

tree_file, out_dir, unaligned_dir, aligned_dir = sys.argv[1:5]
tree = Tree(tree_file)
fasta = tree_file.replace("_alignment.fasta.treefile", "_analysis_candidates.fasta")
seqs = {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}

# Create output directories
os.makedirs(out_dir, exist_ok=True)
os.makedirs(unaligned_dir, exist_ok=True)
os.makedirs(aligned_dir, exist_ok=True)

# Identify clades with support > 70
clade_id = 0
for node in tree.traverse():
    if node.support > 70 and not node.is_leaf():
        clade_seqs = [seqs[leaf.name] for leaf in node.get_leaves() if leaf.name in seqs]
        if clade_seqs:
            SeqIO.write(clade_seqs, f"{out_dir}/clade_{clade_id}.fasta", "fasta")
            SeqIO.write(clade_seqs, f"{unaligned_dir}/clade_{clade_id}.fasta", "fasta")
            clade_id += 1
