#!/usr/bin/env python3
# Purpose: Extract initial sequence-based features from candidate ILPs for ML prediction
# Inputs: analysis/all_candidates.fasta
# Output: analysis/features_initial.csv
# Notes: Structural features deferred to post-ColabFold prediction
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from Bio import SeqIO
import glob
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import yaml

fasta = sys.argv[1]
output_file = sys.argv[2]
with open("config.yaml") as f:
    config = yaml.safe_load(f)

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
    for chunk in pd.read_csv(f, sep="\t", header=None, usecols=[0, 2], chunksize=10000):
        for _, row in chunk.iterrows():
            blast_data[f"{tax_id}_{row[0]}"] = row[2]
features["blast_identity"] = features["id"].map(blast_data).fillna(0)

# Physicochemical properties
features["hydrophobicity"] = [ProteinAnalysis(str(s.seq)).gravy() for s in seqs]
features["charge"] = [ProteinAnalysis(str(s.seq)).charge_at_pH(7.0) for s in seqs]

# InterPro domains
interpro_data = {}
for f in glob.glob("candidates/*_interpro.tsv"):
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", usecols=[0, 4], chunksize=10000):
        for _, row in chunk.groupby(0):
            interpro_data[f"{tax_id}_{row[0].iloc[0]}"] = ";".join(row[4].dropna())
features["domains"] = features["id"].map(interpro_data).fillna("")

features["motifs"] = ""
features.to_csv(output_file, index=False)
