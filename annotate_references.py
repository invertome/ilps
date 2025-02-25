#!/usr/bin/env python3
# Purpose: Annotate UniProt sequences as ILP or non-ILP based on keywords and InterPro domains
# Inputs: 
#   - input/ref_combined.fasta (combined ILP and non-ILP sequences)
#   - input/ref_ILPs_temp_interpro.tsv (InterProScan output)
# Output: input/ref_ILPs.fasta (annotated sequences with [ILP] or [non-ILP] tags)
# Notes: Uses a comprehensive set of keywords and domains for robust annotation
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import pandas as pd
import sys

input_fasta, interpro_tsv, output_fasta = sys.argv[1:4]
records = []

# Load InterPro annotations into a dictionary
interpro_data = pd.read_csv(interpro_tsv, sep="\t", usecols=[0, 4])
interpro_dict = {row[0]: row[4] for _, row in interpro_data.iterrows() if pd.notna(row[4])}

# Define ILP identifiers
ilp_keywords = ["insulin", "igf", "relaxin", "gonadulin", "bombyxin", "insulin-like peptide", "ilp", "neuroparsin"]
ilp_domains = ["PF00049", "PF00159", "PF05389", "PF13715", "IPR022353"]  # Insulin/relaxin family domains

# Annotate each sequence
for rec in SeqIO.parse(input_fasta, "fasta"):
    desc = rec.description.lower()
    interpro = interpro_dict.get(rec.id, "")
    # Check if sequence matches ILP criteria
    is_ilp = any(keyword in desc for keyword in ilp_keywords) or any(domain in interpro for domain in ilp_domains)
    rec.description += " [ILP]" if is_ilp else " [non-ILP]"
    records.append(rec)

# Write annotated sequences to output
SeqIO.write(records, output_fasta, "fasta")
