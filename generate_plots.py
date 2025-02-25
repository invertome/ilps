#!/usr/bin/env python3
# Purpose: Generate multipanel manuscript plots for a given sequence type (prepro, pro, mature)
# Inputs: 
#   - analysis/features.csv (candidate features)
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
#   - clades_ete_dir/ (ETE clade data)
#   - clades_autophy_dir/ (Autophy clade data)
#   - output_dir/ (output directory for plots)
# Outputs: 
#   - output_dir/multipanel_counts_distribution.png (bar and pie charts)
#   - output_dir/multipanel_heatmap_violin.png (heatmap and violin plots)
#   - output_dir/motif_logo_*.png (sequence logos for motifs)
# Notes: Processes data for each sequence type independently; uses taxonkit for taxonomy; chunked processing for memory efficiency
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import subprocess
from matplotlib.gridspec import GridSpec
import glob

# Parse command-line arguments
features_file, pred_file, novel_file, ete_dir, autophy_dir, out_dir = sys.argv[1:7]
features = pd.concat([chunk for chunk in pd.read_csv(features_file, chunksize=10000)])
preds = pd.read_csv(pred_file, names=["ML_Prob"])
novel = pd.read_csv(novel_file, names=["Novel"])

# Add taxonomic information using taxonkit
features["TaxID"] = features["id"].str.split("_").str[0]
tax_cache = f"{out_dir}/tax_map.pkl"
if not os.path.exists(tax_cache):
    tax_map = pd.DataFrame({
        "TaxID": features["TaxID"].unique(),
        "Lineage": [subprocess.check_output(f"taxonkit lineage {taxid}", shell=True).decode().strip() for taxid in features["TaxID"].unique()]
    })
    tax_map[["Phylum", "Class", "Order"]] = tax_map["Lineage"].str.split(";", expand=True)[[1, 2, 3]]
    tax_map.to_pickle(tax_cache)
else:
    tax_map = pd.read_pickle(tax_cache)
features = features.merge(tax_map, on="TaxID")

# Multipanel plot: Counts and distribution by phylum
fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, height_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0, :])
counts = features.groupby(["Phylum", "Novel"]).size().reset_index(name="Candidates")
sns.barplot(x="Phylum", y="Candidates", hue="Novel", data=counts, ax=ax1)
ax1.set_title("ILP Candidates by Phylum")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
dist = features.groupby(["Phylum", "Novel"]).size().unstack(fill_value=0)
dist[0].plot(kind="pie", autopct="%1.1f%%", legend=False, ax=ax2, title="Known ILPs by Phylum")
dist[1].plot(kind="pie", autopct="%1.1f%%", legend=False, ax=ax3, title="Novel ILPs by Phylum")
plt.tight_layout()
os.makedirs(out_dir, exist_ok=True)
plt.savefig(f"{out_dir}/multipanel_counts_distribution.png", bbox_inches="tight")
plt.close()

# Multipanel plot: Heatmap and violin plot
fig = plt.figure(figsize=(14, 12))
gs = GridSpec(2, 1, height_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0])
sns.heatmap(features[["hhblits_prob", "hmm_score", "hhsearch_prob", "tm_score_dynamic_prepro"]], annot=True, cmap="viridis", ax=ax1)
ax1.set_title("Feature Heatmap Across Candidates")
ax2 = fig.add_subplot(gs[1])
sns.violinplot(x="Order", y="hhblits_prob", hue="Novel", data=features, ax=ax2)
ax2.set_title("HHblits Probability Distribution by Order")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{out_dir}/multipanel_heatmap_violin.png", bbox_inches="tight")
plt.close()

# Generate sequence logos for motifs from MEME suite
for meme_file in glob.glob(f"{ete_dir}/meme_*/meme.txt"):
    clade = os.path.basename(os.path.dirname(meme_file)).replace("meme_", "")
    subprocess.run(["ceqlogo", "-i", meme_file, "-o", f"{out_dir}/motif_logo_{clade}.png"])
