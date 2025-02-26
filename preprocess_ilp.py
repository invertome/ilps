#!/usr/bin/env python3
# Purpose: Preprocess ILP sequences into prepro, pro, and mature forms using DeepNeuropePred
# Inputs: input_fasta, sequence_id, type (prepro/pro/mature), output_fasta
# Outputs: output_fasta with the specified form
# Notes: Uses DeepNeuropePred for signal peptide and cleavage site prediction
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import sys
import subprocess
import json
import os
from Bio import SeqIO
import yaml

# Parse command-line arguments
input_fasta, seq_id, seq_type, output_fasta = sys.argv[1:5]

# Load configuration
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Run DeepNeuropePred
output_json = "temp_dnp.json"
deepneuropepred_cmd = [
    "python",
    config["deepneuropepred_script"],
    "--model-file",
    config["deepneuropepred_model"],
    "--input-fasta",
    input_fasta,
    "--output-json",
    output_json,
]
try:
    subprocess.run(deepneuropepred_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error: DeepNeuropePred failed with exit code {e.returncode}", file=sys.stderr)
    sys.exit(1)

# Parse the JSON output
try:
    with open(output_json) as f:
        dnp_predictions = json.load(f)
except (FileNotFoundError, json.JSONDecodeError) as e:
    print(f"Error: Failed to parse DeepNeuropePred output ({e})", file=sys.stderr)
    sys.exit(1)

# Check if the sequence ID exists in the predictions
if seq_id not in dnp_predictions:
    print(f"Error: Sequence ID {seq_id} not found in DeepNeuropePred output", file=sys.stderr)
    sys.exit(1)

prediction = dnp_predictions[seq_id]
signal_pos = prediction["signal_pos"]
cleavage_predictions = prediction["predict"]

# Extract the sequence from the FASTA file
try:
    seq_record = next(rec for rec in SeqIO.parse(input_fasta, "fasta") if rec.id == seq_id)
    full_seq = str(seq_record.seq)
except StopIteration:
    print(f"Error: Sequence ID {seq_id} not found in {input_fasta}", file=sys.stderr)
    sys.exit(1)

# Determine cleavage sites (select the first high-probability site after signal peptide)
cleavage_sites = sorted(
    [pos for pos, prob in cleavage_predictions if prob > 0.5 and pos > signal_pos],
    key=lambda x: x,
)
pro_cleavage = cleavage_sites[0] if cleavage_sites else len(full_seq)  # Default to end if no suitable site

# Split the sequence based on the requested type
if seq_type == "prepro":
    output_seq = full_seq  # Full sequence
elif seq_type == "pro":
    output_seq = full_seq[signal_pos:pro_cleavage]  # From signal peptide end to propeptide cleavage
elif seq_type == "mature":
    output_seq = full_seq[pro_cleavage:]  # From propeptide cleavage to end
else:
    print(f"Error: Invalid sequence type {seq_type}", file=sys.stderr)
    sys.exit(1)

# Write the output to a FASTA file
with open(output_fasta, "w") as f:
    f.write(f">{seq_id}_{seq_type}\n{output_seq}\n")

# Clean up temporary file
os.remove(output_json)
