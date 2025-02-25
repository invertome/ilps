#!/usr/bin/env python3
# Purpose: Combine sequence hits from HHblits and HMMER to identify ILP candidates
# Inputs: 
#   - HHblits output (e.g., candidates/*_hhblits.out)
#   - HMMER output (e.g., candidates/*_hmm.out)
#   - Output FASTA path
# Output: candidates/*_candidates.fasta or preprocess/ref_candidates_temp.fasta (combined hits)
# Notes: Streams sequences to minimize memory usage; takes union of hits
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import sys

hhblits_out, hmm_out, output_fasta = sys.argv[1:4]
hits = set()

# Parse HHblits output for hit IDs
with open(hhblits_out) as fh:
    for line in fh:
        if line.startswith(">"):
            hits.add(line[1:].strip())

# Parse HMMER output for hit IDs
with open(hmm_out) as fh:
    for line in fh:
        if line.startswith(">"):
            hits.add(line[1:].strip())

# Extract sequences from clustered FASTA that match hits
clust_fasta = hhblits_out.replace("_hhblits.out", "_clust_rep_seq.fasta")
with open(output_fasta, "w") as out_fh:
    for rec in SeqIO.parse(clust_fasta, "fasta"):
        if rec.id in hits:
            SeqIO.write(rec, out_fh, "fasta")
