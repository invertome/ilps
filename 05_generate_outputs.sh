#!/bin/bash
# Purpose: Generate tables and plots for manuscript, separated by sequence type (prepro, pro, mature)
# Inputs: 
#   - analysis/features.csv (candidate features)
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
#   - candidates/*_blast.out, candidates/*_interpro.tsv (metadata)
#   - clades_ete_*/, clades_autophy_* (clade data)
# Outputs: 
#   - output/*/*.csv (tables: overview, details, motif enrichment)
#   - output/*/*.png (plots: counts, heatmap, violin, logos)
# Log: pipeline.log
# Notes: Processes prepro, pro, and mature types independently
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f output/.done ]; then
    echo "$(date) - Skipping 05_generate_outputs.sh (already done)" >> pipeline.log
    exit 0
fi

# Set up output directories for each type
mkdir -p output/prepro output/pro output/mature
echo "$(date) - Starting 05_generate_outputs.sh" >> pipeline.log

# Generate tables and plots for each sequence type
python generate_tables.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_prepro/ clades_autophy_prepro/ output/prepro || { echo "$(date) - generate_tables.py failed for prepro" >> pipeline.log; exit 1; }
python generate_tables.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_pro/ clades_autophy_pro/ output/pro || { echo "$(date) - generate_tables.py failed for pro" >> pipeline.log; exit 1; }
python generate_tables.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_mature/ clades_autophy_mature/ output/mature || { echo "$(date) - generate_tables.py failed for mature" >> pipeline.log; exit 1; }

python generate_plots.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_prepro/ clades_autophy_prepro/ output/prepro || { echo "$(date) - generate_plots.py failed for prepro" >> pipeline.log; exit 1; }
python generate_plots.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_pro/ clades_autophy_pro/ output/pro || { echo "$(date) - generate_plots.py failed for pro" >> pipeline.log; exit 1; }
python generate_plots.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_mature/ clades_autophy_mature/ output/mature || { echo "$(date) - generate_plots.py failed for mature" >> pipeline.log; exit 1; }

echo "$(date) - 05_generate_outputs.sh completed" >> pipeline.log
touch output/.done
