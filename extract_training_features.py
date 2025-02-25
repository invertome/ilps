#!/usr/bin/env python3
# Purpose: Extract features from reference ILPs and non-ILPs for ML training
# Inputs: 
#   - preprocess/ref_candidates.fasta (deduplicated candidates)
#   - preprocess/ref_hhblits.out, ref_hmm.out, ref_hhsearch.out, ref_blast.out, ref_interpro.tsv (search outputs)
#   - preprocess/ (directory with PDB files for all types)
# Output: preprocess/ref_features.csv (feature matrix)
# Notes: Includes TM-scores against multiple references for prepro, pro, and mature forms; uses chunked processing for memory efficiency
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from Bio import SeqIO
import glob
import subprocess
import sys
import os
import statistics

fasta, hhblits_out, hmm_out, hhsearch_out, blast_out, interpro_tsv, pdb_dir = sys.argv[1:8]
features = pd.DataFrame()
seqs = list(SeqIO.parse(fasta, "fasta"))
features["id"] = [s.id for s in seqs]

# HHblits probabilities
hhblits_data = {}
with open(hhblits_out) as fh:
    for line in fh:
        if line.startswith(">"):
            id_ = line[1:].strip()
            prob = float(next(fh).split()[2])  # Probability from HHblits
            hhblits_data[id_] = prob
features["hhblits_prob"] = features["id"].map(hhblits_data).fillna(0)

# HMMER scores
hmm_data = {}
with open(hmm_out) as fh:
    for line in fh:
        if line.startswith(">"):
            id_ = line[1:].strip()
            score = float(next(fh).split()[1])  # HMMER score
            hmm_data[id_] = score
features["hmm_score"] = features["id"].map(hmm_data).fillna(0)

# HHsearch probabilities
hhsearch_data = {}
with open(hhsearch_out) as fh:
    for line in fh:
        if "Probab=" in line:
            id_ = line.split()[0]
            prob = float(line.split("Probab=")[1].split()[0])  # HHsearch probability
            hhsearch_data[id_] = prob
features["hhsearch_prob"] = features["id"].map(hhsearch_data).fillna(0)

# BLAST identity (chunked reading for memory efficiency)
blast_chunks = pd.read_csv(blast_out, sep="\t", header=None, usecols=[0, 2], chunksize=10000)
blast_data = {}
for chunk in blast_chunks:
    for _, row in chunk.iterrows():
        blast_data[row[0]] = row[2]  # Percent identity
features["blast_identity"] = features["id"].map(blast_data).fillna(0)

# TM-scores with multiple references for all types
ref_pdbs = {
    "dynamic": "preprocess/ref_structure.pdb",
    "1TRZ": "references/1TRZ.pdb",      # Human insulin
    "6RLX": "references/6RLX.pdb",      # Human relaxin-3
    "bombyxin": "references/bombyxin.pdb"  # Bombyxin-II
}

# Ensure standard reference PDBs are available
os.makedirs("references", exist_ok=True)
if not os.path.exists("references/1TRZ.pdb"):
    subprocess.run(["wget", "-O", "references/1TRZ.pdb", "https://files.rcsb.org/download/1TRZ.pdb"])
if not os.path.exists("references/6RLX.pdb"):
    subprocess.run(["wget", "-O", "references/6RLX.pdb", "https://files.rcsb.org/download/6RLX.pdb"])
if not os.path.exists("references/bombyxin.pdb"):
    with open("references/bombyxin.fasta", "w") as f:
        f.write(">Bombyxin-II\nGIVDECCLRPCSVAALTLREAVNS\n:CLQEGACSVSFGLDAFD")
    subprocess.run(["colabfold_batch", "references/bombyxin.fasta", "references/", "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"], check=True)
    os.rename("references/predict_0.pdb", "references/bombyxin.pdb")

# Generate dynamic reference if missing
if not os.path.exists(ref_pdbs["dynamic"]):
    ilp_seqs = [rec for rec in SeqIO.parse(fasta, "fasta") if "[ILP]" in rec.description]
    median_len = statistics.median([len(rec.seq) for rec in ilp_seqs])
    ref_seq = min(ilp_seqs, key=lambda x: abs(len(x.seq) - median_len))
    with open("preprocess/ref_temp.fasta", "w") as f:
        SeqIO.write(ref_seq, f, "fasta")
    subprocess.run(["colabfold_batch", "preprocess/ref_temp.fasta", "preprocess/", "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"], check=True)
    os.rename("preprocess/predict_0.pdb", ref_pdbs["dynamic"])

# Calculate TM-scores for each structure against all references for all types
tm_data = {ref: {} for ref in ref_pdbs}
for type in ["prepro", "pro", "mature"]:
    for pdb in glob.glob(f"{pdb_dir}/{type}_*.pdb"):
        id_ = os.path.basename(pdb).replace(".pdb", "").replace(f"{type}_", "")
        for ref_name, ref_path in ref_pdbs.items():
            result = subprocess.check_output(["tmalign", ref_path, pdb], text=True)
            tm_score = float([line for line in result.splitlines() if "TM-score=" in line][0].split()[1])
            tm_data[ref_name][id_] = tm_score
    for ref_name in ref_pdbs:
        features[f"tm_score_{ref_name}_{type}"] = features["id"].map(tm_data[ref_name]).fillna(0)

# Add InterPro domains as features (chunked reading)
interpro_chunks = pd.read_csv(interpro_tsv, sep="\t", usecols=[0, 4], chunksize=10000)
interpro_data = {}
for chunk in interpro_chunks:
    for _, row in chunk.groupby(0):
        interpro_data[row[0].iloc[0]] = ";".join(row[4].dropna())
features["domains"] = features["id"].map(interpro_data).fillna("")

# Placeholder for motifs (filled later in analysis)
features["motifs"] = ""
# Write feature matrix to CSV
features.to_csv(sys.argv[8], index=False)
