#!/bin/bash
# Purpose: Generate tables, plots, FASTA files, metadata TSV, and 3D figures for manuscript
# Inputs: analysis/features.csv, analysis/predictions.csv, analysis/novel_candidates.csv, candidates/*_blast.out, candidates/*_interpro.tsv, clades_*
# Outputs: output/*/*.csv, output/*/*.png, output/*/ilps.fasta, output/comparative_metadata.tsv, output/figures/*.png
# Config: config.yaml (max_cpus)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Added PyMOL dependency check
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 05_generate_outputs.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f output/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 05_generate_outputs.sh (already done)" >> pipeline.log
    exit 0
fi

# Check PyMOL dependency
if ! command -v pymol >/dev/null 2>&1; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: PyMOL not installed" >> pipeline.log
    exit 1
fi

mkdir -p output/prepro output/pro output/mature output/figures

for type in prepro pro mature; do
    python generate_tables.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_${type}/ clades_autophy_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_tables.py failed for ${type}" >> pipeline.log; exit 1; }
    python generate_plots.py analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_${type}/ clades_autophy_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_plots.py failed for ${type}" >> pipeline.log; exit 1; }
done

python generate_output_fasta_and_metadata.py analysis/ilp_candidates.fasta analysis/ref_ILPs_filtered.fasta analysis/features.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_interpro.tsv input/[0-9]*_*.fasta output/ || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_output_fasta_and_metadata.py failed" >> pipeline.log; exit 1; }
python generate_structure_figures.py output/figures/ preprocess/ analysis/pdbs/ || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_structure_figures.py failed" >> pipeline.log; exit 1; }

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 05_generate_outputs.sh completed in ${runtime}s" >> pipeline.log
touch output/.done
