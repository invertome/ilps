#!/usr/bin/env python3
# Purpose: Extract sequence-based features from reference ILPs and non-ILPs for ML training
# Inputs: preprocess/ref_candidates.fasta, search outputs
# Output: preprocess/ref_features.csv
# Notes: Structural features deferred to candidate prediction
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
from Bio import SeqIO
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import yaml

fasta, hhblits_out, hmm_out, hhsearch_out, blast_out, interpro_tsv, _ = sys.argv[1:8]  # Last arg unused
with open("config.yaml") as f:
    config = yaml.safe_load(f)

features = pd.DataFrame()
seqs = list(SeqIO.parse(fasta, "fasta"))
features["id"] = [s.id for s in seqs]

# HHblits probabilities
hhblits_data = {}
with open(hhblits_out) as fh:
    for line in fh:
        if line.startswith(">"):
            id_ = line[1:].strip()
            prob = float(next(fh).split()[2])
            hhblits_data[id_] = prob
features["hhblits_prob"] = features["id"].map(hhblits_data).fillna(0)

# HMMER scores
hmm_data = {}
with open(hmm_out) as fh:
    for line in fh:
        if line.startswith(">"):
            id_ = line[1:].strip()
            score = float(next(fh).split()[1])
            hmm_data[id_] = score
features["hmm_score"] = features["id"].map(hmm_data).fillna(0)

# HHsearch probabilities
hhsearch_data = {}
with open(hhsearch_out) as fh:
    for line in fh:
        if "Probab=" in line:
            id_ = line.split()[0]
            prob = float(line.split("Probab=")[1].split()[0])
            hhsearch_data[id_] = prob
features["hhsearch_prob"] = features["id"].map(hhsearch_data).fillna(0)

# BLAST identity
blast_chunks = pd.read_csv(blast_out, sep="\t", header=None, usecols=[0, 2], chunksize=10000)
blast_data = {}
for chunk in blast_chunks:
    for _, row in chunk.iterrows():
        blast_data[row[0]] = row[2]
features["blast_identity"] = features["id"].map(blast_data).fillna(0)

# Physicochemical properties
features["hydrophobicity"] = [ProteinAnalysis(str(s.seq)).gravy() for s in seqs]
features["charge"] = [ProteinAnalysis(str(s.seq)).charge_at_pH(7.0) for s in seqs]

# InterPro domains
interpro_chunks = pd.read_csv(interpro_tsv, sep="\t", usecols=[0, 4], chunksize=10000)
interpro_data = {}
for chunk in interpro_chunks:
    for _, row in chunk.groupby(0):
        interpro_data[row[0].iloc[0]] = ";".join(row[4].dropna())
features["domains"] = features["id"].map(interpro_data).fillna("")

features["motifs"] = ""
features.to_csv(sys.argv[8], index=False)
