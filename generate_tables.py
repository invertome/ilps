#!/usr/bin/env python3
# Purpose: Generate manuscript tables for a given sequence type (prepro, pro, mature)
# Inputs: 
#   - analysis/features.csv (candidate features)
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
#   - candidates/*_blast.out (BLAST results)
#   - candidates/*_interpro.tsv (InterPro annotations)
#   - clades_ete_dir/ (ETE clade data)
#   - clades_autophy_dir/ (Autophy clade data)
#   - output_dir/ (output directory for tables)
# Outputs: 
#   - output_dir/ilp_overview.csv (summary by transcriptome)
#   - output_dir/ilp_details.csv (per-sequence details)
#   - output_dir/motif_enrichment.csv (motif enrichment analysis)
# Notes: Processes data for each sequence type independently; uses chunked processing for memory efficiency
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
import glob
import os
import sys

# Parse command-line arguments
features_file, pred_file, novel_file, *extra_files, out_dir = sys.argv[1:-1]
features = pd.concat([chunk for chunk in pd.read_csv(features_file, chunksize=10000)])
preds = pd.read_csv(pred_file, names=["ML_Prob"])
novel = pd.read_csv(novel_file, names=["Novel"])
blast_files = [f for f in extra_files if f.endswith("_blast.out")]
interpro_files = [f for f in extra_files if f.endswith("_interpro.tsv")]
ete_dir, autophy_dir = extra_files[-2:]

# Overview table: summarize candidates by transcriptome
overview = features.groupby(features["id"].str.split("_").str[1]).agg({
    "id": "count",
    "ML_Prob": "mean",
    "hhblits_prob": "mean",
    "tm_score_dynamic_prepro": "mean",
    "tm_score_1TRZ_prepro": "mean",
    "tm_score_6RLX_prepro": "mean",
    "tm_score_bombyxin_prepro": "mean"
}).rename(columns={
    "id": "Candidates",
    "ML_Prob": "Avg_ML_Prob",
    "hhblits_prob": "Avg_HHblits_Prob",
    "tm_score_dynamic_prepro": "Avg_TM_Score_Dynamic",
    "tm_score_1TRZ_prepro": "Avg_TM_Score_1TRZ",
    "tm_score_6RLX_prepro": "Avg_TM_Score_6RLX",
    "tm_score_bombyxin_prepro": "Avg_TM_Score_Bombyxin"
})
overview["Known_ILPs"] = ((preds["ML_Prob"] > 0.7) & (~novel["Novel"].astype(bool))).groupby(features["id"].str.split("_").str[1]).sum()
overview["Novel_ILPs"] = novel["Novel"].groupby(features["id"].str.split("_").str[1]).sum()
overview = overview.reset_index(names="Transcriptome")

# Details table: per-sequence annotations
blast_hits = {}
for f in blast_files:
    tax_id = os.path.basename(f).split("_")[0]
    for chunk in pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], chunksize=10000):
        for _, row in chunk.iterrows():
            blast_hits[f"{tax_id}_{row[0]}"] = f"{row[1]} ({row[2]}%)"
details = pd.concat([features, preds, novel], axis=1)
details["Top_BLAST_Hit"] = details["id"].map(blast_hits).fillna("No Hit")
details["Clade"] = details["motifs"].apply(lambda x: x.split(";")[0] if x else "Unclustered")

# Motif enrichment table: analyze motifs from MEME suite
motif_data = []
for clade_dir in [ete_dir, autophy_dir]:
    for ame_file in glob.glob(f"{clade_dir}/ame_*.txt"):
        clade = os.path.basename(ame_file).replace("ame_", "").replace(".txt", "")
        df = pd.read_csv(ame_file, sep="\t", comment="#", usecols=["motif_id", "p-value"])
        fimo_file = ame_file.replace("ame_", "fimo_")
        fimo_hits = len(pd.read_csv(fimo_file, sep="\t", comment="#")) if os.path.exists(fimo_file) else 0
        for _, row in df.iterrows():
            motif_data.append({
                "Clade": clade,
                "Motif": row["motif_id"],
                "P-value": row["p-value"],
                "FIMO_Hits": fimo_hits
            })
motif_table = pd.DataFrame(motif_data)

# Ensure output directory exists and write tables
os.makedirs(out_dir, exist_ok=True)
overview.to_csv(f"{out_dir}/ilp_overview.csv", index=False)
details.to_csv(f"{out_dir}/ilp_details.csv", index=False)
motif_table.to_csv(f"{out_dir}/motif_enrichment.csv", index=False)
