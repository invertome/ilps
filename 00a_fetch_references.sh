#!/bin/bash
# Purpose: Fetch ILP and non-ILP reference sequences from UniProt and curated sources for a balanced training set
# Inputs: None (queries UniProt API and curated FASTA)
# Outputs: input/ref_ILPs.fasta (annotated ILP/non-ILP sequences)
# Config: config.yaml (max_cpus, thresholds)
# Log: pipeline.log (tracks progress, errors, and profiling)
# Notes: Ensures diverse ILP representation with high-quality curated references
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

# Load config
max_cpus=$(yq e '.max_cpus' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 00a_fetch_references.sh" >> pipeline.log

# Profiling setup with psutil
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f input/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 00a_fetch_references.sh (already done)" >> pipeline.log
    exit 0
fi

# Create input directory
mkdir -p input

# Fetch ILP sequences with broad keywords across metazoans
curl -o input/ref_ILPs_raw.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(reviewed:true)%20AND%20(taxonomy:%22Metazoa%20[33208]%22)&format=fasta" \
    || { echo "$(date '+%Y-%m-%d %H:%M:%S') - UniProt ILP query failed" >> pipeline.log; exit 1; }

# Fetch curated high-quality ILP references (e.g., from literature, Swiss-Prot)
curl -o input/ref_ILPs_curated.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(accession:P01308%20OR%20accession:P01315%20OR%20accession:P01344%20OR%20accession:P51435)%20AND%20(reviewed:true)&format=fasta" \
    || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Uni
