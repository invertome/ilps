#!/usr/bin/env python3
# Purpose: Extract features from candidate ILPs for ML prediction
# Inputs: analysis/all_candidates.fasta (combined candidates)
# Output: analysis/features.csv (feature matrix)
# Notes: Similar to training features but uses candidate-specific outputs; includes TM-scores for all types; uses chunked processing
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from Bio import SeqIO
import glob
import subprocess
import sys
import os
import statistics

fasta = sys.argv[1]
output_file = sys.argv[2]
features = pd.DataFrame()
seqs = list(SeqIO.parse(fasta, "fasta"))
features["id"] = [s.id for s in seqs]

# HHblits probabilities from candidate searches
hhblits_data = {}
for f in glob.glob("candidates/*_hhblits.out"):
    tax_id = os.path.basename(f).split("_")[0]
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"):
                id_ = line[1:].strip()
                prob = float(next(fh).split()[2])
                hhblits_data[f"{tax_id}_{id_}"] = prob
features["hhblits_prob"] = features["id"].map(hhblits_data).fillna(0)

# HMMER scores from candidate searches
hmm_data = {}
for f in glob.glob("candidates/*_hmm.out"):
    tax_id = os.path.basename(f).split("_")[0]
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"):
                id_ = line[1:].strip()
                score = float(next(fh).split()[1])
                hmm_data[f"{tax_id}_{id_}"] = score
features["hmm_score"] = features["id"].map(hmm_data).fillna(0)

# HHsearch probabilities from candidate searches
hhsearch_data = {}
for f in glob.glob("candidates/*_hhsearch.out"):
    tax_id = os.path.basename(f).split("_")[0]
    with open(f) as fh:
        for line in fh:
            if "Probab=" in line:
                id_ = line.split()[0]
                prob = float(line.split("Probab=")[1].split()[0])
                hhsearch_data[f"{tax_id}_{id_}"] = prob
features["hhsearch_prob"] = features["id"].map(hhsearch_data).fillna(0)

# BLAST identity from candidate searches (chunked reading)
blast_data = {}
for f in glob.glob("candidates/*_blast.out"):
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", header=None, usecols=[0, 2], chunksize=10000):
        for _, row in chunk.iterrows():
            blast_data[f"{tax_id}_{row[0]}"] = row[2]
features["blast_identity"] = features["id"].map(blast_data).fillna(0)

# TM-scores with multiple references for all types
ref_pdbs = {
    "dynamic": "preprocess/ref_structure.pdb",
    "1TRZ": "references/1TRZ.pdb",
    "6RLX": "references/6RLX.pdb",
    "bombyxin": "references/bombyxin.pdb"
}

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

if not os.path.exists(ref_pdbs["dynamic"]) and os.path.exists("analysis/ref_ILPs_filtered.fasta"):
    ilp_seqs = list(SeqIO.parse("analysis/ref_ILPs_filtered.fasta", "fasta"))
    median_len = statistics.median([len(rec.seq) for rec in ilp_seqs])
    ref_seq = min(ilp_seqs, key=lambda x: abs(len(x.seq) - median_len))
    with open("preprocess/ref_temp.fasta", "w") as f:
        SeqIO.write(ref_seq, f, "fasta")
    subprocess.run(["colabfold_batch", "preprocess/ref_temp.fasta", "preprocess/", "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"], check=True)
    os.rename("preprocess/predict_0.pdb", ref_pdbs["dynamic"])

tm_data = {ref: {} for ref in ref_pdbs}
for type in ["prepro", "pro", "mature"]:
    for pdb in glob.glob(f"analysis/pdbs/{type}_*.pdb"):
        id_ = os.path.basename(pdb).replace(".pdb", "").replace(f"{type}_", "")
        for ref_name, ref_path in ref_pdbs.items():
            result = subprocess.check_output(["tmalign", ref_path, pdb], text=True)
            tm_score = float([line for line in result.splitlines() if "TM-score=" in line][0].split()[1])
            tm_data[ref_name][id_] = tm_score
    for ref_name in ref_pdbs:
        features[f"tm_score_{ref_name}_{type}"] = features["id"].map(tm_data[ref_name]).fillna(0)

# Add InterPro domains as features (chunked reading)
interpro_data = {}
for f in glob.glob("candidates/*_interpro.tsv"):
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", usecols=[0, 4], chunksize=10000):
        for _, row in chunk.groupby(0):
            interpro_data[f"{tax_id}_{row[0].iloc[0]}"] = ";".join(row[4].dropna())
features["domains"] = features["id"].map(interpro_data).fillna("")

# Placeholder for motifs (filled later)
features["motifs"] = ""
features.to_csv(output_file, index=False)
