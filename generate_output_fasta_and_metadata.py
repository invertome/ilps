#!/usr/bin/env python3
# Purpose: Generate FASTA files and metadata TSV for predicted ILPs and references
# Inputs: ilp_candidates.fasta, ref_ILPs_filtered.fasta, features.csv, predictions.csv, novel_candidates.csv, interpro files, input FASTA files, output dir
# Outputs: output/*/ilps.fasta, output/comparative_metadata.tsv
# Notes: Includes detailed headers and comprehensive metadata
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

from Bio import SeqIO
import pandas as pd
import glob
import sys
import os
import subprocess

ilp_fasta, ref_fasta, features_file, pred_file, novel_file, *extra_files, output_dir = sys.argv[1:]
interpro_files = [f for f in extra_files if f.endswith("_interpro.tsv")]
input_files = [f for f in extra_files if f.startswith("input/")]

features = pd.read_csv(features_file)
preds = pd.read_csv(pred_file, names=["ML_Prob"])
novel = pd.read_csv(novel_file, names=["Novel"])
ilp_seqs = list(SeqIO.parse(ilp_fasta, "fasta"))
ref_seqs = list(SeqIO.parse(ref_fasta, "fasta"))

# InterPro domains
interpro_data = {}
for f in interpro_files:
    tax_id = os.path.basename(f).split("_")[0]
    df = pd.read_csv(f, sep="\t", usecols=[0, 4])
    for _, row in df.groupby(0):
        interpro_data[f"{tax_id}_{row[0].iloc[0]}"] = ";".join(row[4].dropna()) if len(row[4].dropna()) > 0 else ""

# Taxonomy and nucleotide sequences
tax_data = {}
nuc_seqs = {}
for f in input_files:
    tax_id = os.path.basename(f).split("_")[0]
    lineage = subprocess.check_output(f"taxonkit lineage {tax_id}", shell=True).decode().strip().split(";")
    tax_data[tax_id] = {rank.split(" ")[0]: rank.split(" ")[1][1:-1] for rank in lineage if len(rank.split(" ")) > 1}
    for rec in SeqIO.parse(f, "fasta"):
        if f.endswith("_preprocessed.fasta"):
            nuc_seqs[rec.id] = str(rec.seq)
        else:
            for aa_rec in ilp_seqs:
                if aa_rec.id.startswith(tax_id):
                    nuc_seqs[aa_rec.id] = str(rec.seq)  # Simplified mapping

# Generate FASTA files
for type in ["prepro", "pro", "mature"]:
    with open(f"{output_dir}/{type}/ilps.fasta", "w") as f:
        for seq, prob, nov in zip(ilp_seqs + ref_seqs, preds["ML_Prob"] + [1] * len(ref_seqs), novel["Novel"].tolist() + [0] * len(ref_seqs)):
            tax_id = seq.id.split("_")[0] if "_" in seq.id else "ref"
            domains = interpro_data.get(seq.id, "")
            header = f">{seq.id}_{type} ML_Prob={prob:.3f} Novel={nov} Domains={domains} Taxonomy={'|'.join([f'{k}={v}' for k, v in tax_data.get(tax_id, {}).items()])}"
            aa_seq = next((str(s.seq) for s in SeqIO.parse(f"analysis/pdbs/{type}_{seq.id}.fasta", "fasta")), str(seq.seq))
            f.write(f"{header}\n{aa_seq}\n")

# Generate metadata TSV
metadata = pd.DataFrame({
    "SequenceID": [s.id for s in ilp_seqs + ref_seqs],
    **{f"Tax_{rank}": [tax_data.get(s.id.split("_")[0] if "_" in s.id else "ref", {}).get(rank, "") for s in ilp_seqs + ref_seqs] for rank in ["phylum", "class", "order", "family", "genus", "species"]},
    "Annotations": ["ILP" if p > 0.7 else "Non-ILP" for p in preds["ML_Prob"]] + ["ILP"] * len(ref_seqs),
    "Clades": [f["motifs"].split(";")[0] if f["motifs"] else "Unclustered" for _, f in features.iterrows()] + ["Reference"] * len(ref_seqs),
    "Domains": [interpro_data.get(s.id, "") for s in ilp_seqs + ref_seqs],
    "NucleotideSeq": [nuc_seqs.get(s.id, "") for s in ilp_seqs + ref_seqs],
    "AminoAcidSeq": [str(s.seq) for s in ilp_seqs + ref_seqs],
    "GenBankNucID": [""] * (len(ilp_seqs) + len(ref_seqs)),
    "GenBankProtID": [""] * (len(ilp_seqs) + len(ref_seqs))
})
metadata.to_csv(f"{output_dir}/comparative_metadata.tsv", sep="\t", index=False)
