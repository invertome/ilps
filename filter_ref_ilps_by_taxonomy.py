#!/usr/bin/env python3
# Purpose: Filter reference ILPs to match the highest common taxonomic rank
# Inputs: 
#   - input/ref_ILPs.fasta (annotated ILP/non-ILP sequences)
#   - TaxID (e.g., 33208 from common_taxonomy.txt)
#   - Output FASTA path (e.g., analysis/ref_ILPs_filtered.fasta)
# Output: analysis/ref_ILPs_filtered.fasta (filtered ILP sequences)
# Notes: Uses taxonkit to check lineage inclusion
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import sys
import subprocess

ref_fasta, taxid, output_fasta = sys.argv[1:4]
records = []
for rec in SeqIO.parse(ref_fasta, "fasta"):
    if "[ILP]" in rec.description:  # Only keep ILPs
        # Extract UniProt TaxID (e.g., from "sp|P01308|INS_HUMAN")
        uni_taxid = rec.id.split("|")[1]
        lineage = subprocess.check_output(f"taxonkit lineage {uni_taxid}", shell=True).decode().strip()
        if f"[{taxid}]" in lineage:  # Check if TaxID is in lineage
            records.append(rec)
SeqIO.write(records, output_fasta, "fasta")
