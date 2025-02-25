#!/usr/bin/env Rscript
# Purpose: Generate a majority-rule consensus tree with clade support values from multiple input phylogenies
# Inputs: Multiple tree files (e.g., sequence and structural trees), output path
# Output: Consensus tree file with support values (frequency of clades) as node labels
# Notes: Uses ape package; majority-rule with p=0.5; calculates clade frequencies
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

library(ape)

# Function to check if a clade exists in a tree
match.clade <- function(tree, clade) {
    tips <- clade$tip.label
    if (length(tips) <= 1) return(TRUE)  # Single tip is always present
    mrca <- getMRCA(tree, tips)
    if (is.null(mrca)) return(FALSE)
    desc <- extract.clade(tree, mrca)$tip.label
    return(setequal(tips, desc))
}

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)
input_trees <- args[1:(length(args)-1)]  # All but last argument are input trees
output_file <- args[length(args)]        # Last argument is output file

# Read input trees
trees <- lapply(input_trees, read.tree)

# Generate majority-rule consensus tree (p=0.5)
cons_tree <- consensus(trees, p=0.5, check.labels=TRUE)

# Add support values (frequency of clades across input trees)
cons_tree$node.label <- sapply(1:(length(cons_tree$tip.label)-1), function(i) {
    clade <- extract.clade(cons_tree, i + length(cons_tree$tip.label))
    mean(sapply(trees, function(t) match.clade(t, clade)))
})

# Write consensus tree to file
write.tree(cons_tree, output_file)
