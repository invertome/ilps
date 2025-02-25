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
    || { echo "$(date '+%Y-%m-%d %H:%M:%S') - UniProt curated ILP query failed" >> pipeline.log; exit 1; }

# Combine curated and raw ILP sequences
cat input/ref_ILPs_raw.fasta input/ref_ILPs_curated.fasta > input/ref_ILPs_temp.fasta

# Annotate with InterProScan
interpro_path=$(yq e '.interpro_path' config.yaml)
"$interpro_path" -i input/ref_ILPs_temp.fasta -dp -f tsv -iprlookup -goterms -pa -cpu "$max_cpus" > input/ref_ILPs_temp_interpro.tsv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - InterProScan failed" >> pipeline.log; exit 1; }

# Count ILPs for balancing
ilp_count=$(grep -c "^>" input/ref_ILPs_temp.fasta)

# Fetch random non-ILP sequences
curl -o input/ref_non_ILPs.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)%20AND%20(taxonomy:%22Metazoa%20[33208]%22)%20AND%20NOT%20(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(length:[50%20TO%20200])&format=fasta&size=${ilp_count}&sort=rand" \
    || { echo "$(date '+%Y-%m-%d %H:%M:%S') - UniProt non-ILP query failed" >> pipeline.log; exit 1; }

# Combine and annotate
cat input/ref_ILPs_temp.fasta input/ref_non_ILPs.fasta > input/ref_combined.fasta
python annotate_references.py input/ref_combined.fasta input/ref_ILPs_temp_interpro.tsv input/ref_ILPs.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - annotate_references.py failed" >> pipeline.log; exit 1; }

# Clean up
pigz -f input/ref_ILPs_raw.fasta input/ref_ILPs_curated.fasta input/ref_non_ILPs.fasta input/ref_combined.fasta input/ref_ILPs_temp.fasta input/ref_ILPs_temp_interpro.tsv

# Profiling and completion
end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 00a_fetch_references.sh completed in ${runtime}s" >> pipeline.log
touch input/.done
