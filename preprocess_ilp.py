#!/usr/bin/env python3
# Purpose: Process ILP sequences into prepropeptides, propeptides, or mature peptides using SignalP and InterPro
# Inputs: 
#   - Input FASTA (e.g., preprocess/ref_candidates.fasta)
#   - Sequence ID
#   - Type (prepro, pro, mature)
#   - Output FASTA path (e.g., preprocess/prepro_${seq}.fasta)
# Output: FASTA file with processed sequence
# Notes: Uses SignalP for signal peptide cleavage and InterPro for mature peptide boundaries; includes error handling
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import subprocess
import sys
import os
import tempfile
import pandas as pd

fasta, seq_id, seq_type, output_fasta = sys.argv[1:5]
record = next(r for r in SeqIO.parse(fasta, "fasta") if r.id == seq_id)
seq = str(record.seq)

# Create temporary FASTA file for SignalP
with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_fasta:
    temp_fasta.write(f">{seq_id}\n{seq}")
    temp_fasta_path = temp_fasta.name

# Run SignalP to predict signal peptide cleavage site with error handling
signalp_out = f"{temp_fasta_path}.signalp"
try:
    subprocess.run(["signalp6", "-fasta", temp_fasta_path, "-org", "euk", "-format", "short", "-prefix", temp_fasta_path], check=True)
except subprocess.CalledProcessError:
    print(f"SignalP failed for {seq_id}", file=sys.stderr)
    os.remove(temp_fasta_path)
    sys.exit(1)

# Parse SignalP output for cleavage position
cleavage_pos = 0
with open(f"{signalp_out}_prediction_results.txt") as f:
    lines = f.readlines()[1:]  # Skip header
    for line in lines:
        if line.startswith(seq_id):
            parts = line.split()
            if parts[1] == "SP":  # Signal peptide predicted
                cleavage_pos = int(parts[2])  # Position after cleavage
            break

# Process sequence based on type
if seq_type == "prepro":
    processed_seq = seq  # Full prepropeptide
elif seq_type == "pro" or seq_type == "mature":
    # Remove signal peptide if predicted
    pro_seq = seq[cleavage_pos:] if cleavage_pos > 0 else seq
    if seq_type == "pro":
        processed_seq = pro_seq  # Propeptide (no further cleavage)
    else:
        # Run InterProScan to identify mature peptide regions (e.g., insulin family domains) with error handling
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_fasta:
            temp_fasta.write(f">{seq_id}\n{pro_seq}")
            temp_fasta_path = temp_fasta.name
        interpro_out = f"{temp_fasta_path}.interpro"
        try:
            subprocess.run(["interproscan.sh", "-i", temp_fasta_path, "-dp", "-f", "tsv", "-iprlookup", "-goterms", "-pa", "-o", interpro_out], check=True)
        except subprocess.CalledProcessError:
            print(f"InterProScan failed for {seq_id}", file=sys.stderr)
            os.remove(temp_fasta_path)
            sys.exit(1)
        
        # Parse InterPro output for insulin family domain (PF00049)
        df = pd.read_csv(interpro_out, sep="\t", usecols=[0, 4, 6, 7], names=["id", "domain", "start", "end"], skiprows=1)
        mature_region = df[df["domain"] == "PF00049"]
        if not mature_region.empty:
            start = mature_region["start"].iloc[0] - 1  # 0-based
            end = mature_region["end"].iloc[0]
            mature_seq = pro_seq[start:end]
            # Split into A- and B-chains based on cysteine pairs
            cysteines = [i for i, c in enumerate(mature_seq) if c == "C"]
            if len(cysteines) >= 4:
                b_chain = mature_seq[:cysteines[1] + 1]  # Up to second cysteine
                a_chain = mature_seq[cysteines[-2]:]      # Last two cysteines onward
                processed_seq = f"{b_chain}:{a_chain}"    # Multi-chain format
            else:
                processed_seq = mature_seq  # Fallback if not enough cysteines
        else:
            processed_seq = pro_seq  # Fallback if no domain found
        
        # Clean up InterPro temp files
        os.remove(temp_fasta_path)
        os.remove(interpro_out)
else:
    raise ValueError("Invalid sequence type. Use 'prepro', 'pro', or 'mature'.")

# Write processed sequence to output FASTA
with open(output_fasta, "w") as f:
    f.write(f">{seq_id}_{seq_type}\n{processed_seq}\n")

# Clean up SignalP temp files
os.remove(temp_fasta_path)
os.remove(f"{signalp_out}_prediction_results.txt")
os.remove(f"{signalp_out}_summary.signalp6")
