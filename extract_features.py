#!/usr/bin/env python3
# Purpose: Extract final features (sequence + structural) from candidate ILPs for ML prediction and output
# Inputs: analysis/all_candidates.fasta, analysis/pdbs/*.pdb
# Output: analysis/features_final.csv
# Notes: Includes structural features post-ColabFold prediction
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from Bio import SeqIO
import glob
import subprocess
import sys
import os
import statistics
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import yaml
import psutil

fasta = sys.argv[1]
output_file = sys.argv[2]
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Dynamic chunk size based on available memory
chunksize = int(psutil.virtual_memory().available / (1024 * 1024 * 10))  # Adjust as needed

features = pd.DataFrame()
seqs = list(SeqIO.parse(fasta, "fasta"))
features["id"] = [s.id for s in seqs]

# HHblits probabilities
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

# HMMER scores
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

# HHsearch probabilities
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

# BLAST identity
blast_data = {}
for f in glob.glob("candidates/*_blast.out"):
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", header=None, usecols=[0, 2], chunksize=chunksize):
        for _, row in chunk.iterrows():
            blast_data[f"{tax_id}_{row[0]}"] = row[2]
features["blast_identity"] = features["id"].map(blast_data).fillna(0)

# Physicochemical properties
features["hydrophobicity"] = [ProteinAnalysis(str(s.seq)).gravy() for s in seqs]
features["charge"] = [ProteinAnalysis(str(s.seq)).charge_at_pH(7.0) for s in seqs]

# TM-scores and pLDDT scores from candidate PDBs
ref_pdbs = config["tmalign_refs"]
os.makedirs("references", exist_ok=True)
for ref_name, ref_path in ref_pdbs.items():
    if not os.path.exists(ref_path) and ref_name != "dynamic":
        if ref_name == "1TRZ":
            subprocess.run(["wget", "-O", ref_path, "https://files.rcsb.org/download/1TRZ.pdb"])
        elif ref_name == "6RLX":
            subprocess.run(["wget", "-O", ref_path, "https://files.rcsb.org/download/6RLX.pdb"])
        elif ref_name == "bombyxin":
            with open("references/bombyxin.fasta", "w") as f:
                f.write(">Bombyxin-II\nGIVDECCLRPCSVAALTLREAVNS\n:CLQEGACSVSFGLDAFD")
            subprocess.run([config["colabfold_path"], "references/bombyxin.fasta", "references/", "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"])
            os.rename("references/predict_0.pdb", ref_path)

if not os.path.exists(ref_pdbs["dynamic"]) and os.path.exists("analysis/ref_ILPs_filtered.fasta"):
    ilp_seqs = list(SeqIO.parse("analysis/ref_ILPs_filtered.fasta", "fasta"))
    median_len = statistics.median([len(rec.seq) for rec in ilp_seqs])
    ref_seq = min(ilp_seqs, key=lambda x: abs(len(x.seq) - median_len))
    with open("preprocess/ref_temp.fasta", "w") as f:
        SeqIO.write(ref_seq, f, "fasta")
    subprocess.run([config["colabfold_path"], "preprocess/ref_temp.fasta", "preprocess/", "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"])
    os.rename("preprocess/predict_0.pdb", ref_pdbs["dynamic"])

tm_data = {ref: {} for ref in ref_pdbs}
plddt_data = {ref: {} for ref in ref_pdbs}
for type in ["prepro", "pro", "mature"]:
    for pdb in glob.glob(f"analysis/pdbs/{type}_*.pdb"):
        id_ = os.path.basename(pdb).replace(".pdb", "").replace(f"{type}_", "")
        for ref_name, ref_path in ref_pdbs.items():
            result = subprocess.check_output(["tmalign", ref_path, pdb], text=True)
            tm_score = float([line for line in result.splitlines() if "TM-score=" in line][0].split()[1])
            tm_data[ref_name][id_] = tm_score
            with open(pdb) as f:
                plddt = statistics.mean([float(line.split()[10]) for line in f if line.startswith("ATOM")])
            plddt_data[ref_name][id_] = plddt
    for ref_name in ref_pdbs:
        features[f"tm_score_{ref_name}_{type}"] = features["id"].map(tm_data[ref_name]).fillna(0)
        features[f"plddt_{ref_name}_{type}"] = features["id"].map(plddt_data[ref_name]).fillna(0)

# InterPro domains
interpro_data = {}
for f in glob.glob("candidates/*_interpro.tsv"):
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", usecols=[0, 4], chunksize=chunksize):
        for _, row in chunk.groupby(0):
            interpro_data[f"{tax_id}_{row[0].iloc[0]}"] = ";".join(row[4].dropna())
features["domains"] = features["id"].map(interpro_data).fillna("")

features["motifs"] = ""
features.to_csv(output_file, index=False)
