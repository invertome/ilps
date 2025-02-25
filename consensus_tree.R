#!/usr/bin/env Rscript
# Purpose: Generate a consensus tree from multiple input phylogenies
# Inputs: Multiple tree files (e.g., analysis/prepro_alignment.fasta.treefile, analysis/pro_alignment.fasta.treefile), output path
# Output: Consensus tree file (e.g., analysis/prepro_pro_seq_consensus.tre)
# Notes: Uses ape package; takes majority rule consensus with p=0.5; handles both sequence and structural trees
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

library(ape)
args <- commandArgs(trailingOnly=TRUE)
input_trees <- args[1:(length(args)-1)]  # All but last argument are input trees
output_file <- args[length(args)]        # Last argument is output file

# Read input trees
trees <- lapply(input_trees, read.tree)

# Generate consensus tree with majority rule (p=0.5)
cons_tree <- consensus(trees, p=0.5, check.labels=TRUE)

# Write consensus tree to file
write.tree(cons_tree, output_file)
