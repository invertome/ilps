#!/usr/bin/env python3
# Purpose: Determine the highest common taxonomic rank among input transcriptomes
# Inputs: 
#   - input/*.fasta (transcriptome files with TaxID prefixes, e.g., 9606_T1.fasta)
#   - Output file path (e.g., analysis/common_taxonomy.txt)
# Output: analysis/common_taxonomy.txt (rank and TaxID, e.g., "class 33208")
# Notes: Uses taxonkit for lineage lookup with caching to improve efficiency
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import sys
import subprocess
import pandas as pd
import os
import pickle

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]
taxids = [f.split("_")[0] for f in input_files]  # Extract TaxIDs from filenames

# Load or create taxonomy cache
cache_file = "taxonomy_cache.pkl"
if os.path.exists(cache_file):
    with open(cache_file, "rb") as f:
        tax_cache = pickle.load(f)
else:
    tax_cache = {}

# Get lineages for each TaxID
lineages = []
for taxid in taxids:
    if taxid not in tax_cache:
        lineage = subprocess.check_output(f"taxonkit lineage {taxid}", shell=True).decode().strip()
        tax_cache[taxid] = lineage
    lineages.append(tax_cache[taxid].split(";"))

# Save updated cache
with open(cache_file, "wb") as f:
    pickle.dump(tax_cache, f)

# Find highest common rank
common_rank = None
common_taxid = None
for i in range(min(len(l) for l in lineages)):
    if all(lineages[0][i] == l[i] for l in lineages[1:]):
        rank_parts = lineages[0][i].split(" ")
        common_rank = rank_parts[0]  # e.g., "class"
        common_taxid = rank_parts[1][1:-1]  # e.g., "33208"
    else:
        break

# Write result to file
with open(output_file, "w") as f:
    f.write(f"{common_rank} {common_taxid}")
