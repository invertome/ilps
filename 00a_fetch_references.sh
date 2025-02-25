#!/bin/bash
# Purpose: Fetch ILP and non-ILP reference sequences from UniProt to build a balanced training set
# Inputs: None (queries UniProt API directly)
# Outputs: input/ref_ILPs.fasta (annotated ILP and non-ILP sequences)
# Log: pipeline.log (tracks progress and errors)
# Notes: Ensures diverse ILP representation and equal non-ILP count for ML training
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f input/.done ]; then
    echo "$(date) - Skipping 00a_fetch_references.sh (already done)" >> pipeline.log
    exit 0
fi

# Create input directory if it doesn't exist
mkdir -p input
echo "$(date) - Starting 00a_fetch_references.sh" >> pipeline.log

# Fetch ILP sequences with broad keywords across metazoans
curl -o input/ref_ILPs_raw.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(reviewed:true)%20AND%20(taxonomy:%22Metazoa%20[33208]%22)&format=fasta" \
    || { echo "$(date) - UniProt ILP query failed" >> pipeline.log; exit 1; }
curl -o input/ref_ILPs_model.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(reviewed:true)%20AND%20((taxonomy:%22Homo%20sapiens%20[9606]%22)%20OR%20(taxonomy:%22Mus%20musculus%20[10090]%22)%20OR%20(taxonomy:%22Drosophila%20melanogaster%20[7227]%22))&format=fasta" \
    || { echo "$(date) - UniProt model ILP query failed" >> pipeline.log; exit 1; }
curl -o input/ref_ILPs_nonmodel.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(reviewed:true)%20AND%20((taxonomy:%22Branchiostoma%20[7739]%22)%20OR%20(taxonomy:%22Strongylocentrotus%20purpuratus%20[7668]%22)%20OR%20(taxonomy:%22Cnidaria%20[6073]%22))&format=fasta" \
    || { echo "$(date) - UniProt non-model ILP query failed" >> pipeline.log; exit 1; }
cat input/ref_ILPs_raw.fasta input/ref_ILPs_model.fasta input/ref_ILPs_nonmodel.fasta > input/ref_ILPs_temp.fasta

# Annotate with InterProScan for domain-based ILP identification
max_cpus=$(nproc)
interproscan.sh -i input/ref_ILPs_temp.fasta -dp -f tsv -iprlookup -goterms -pa -cpu "$max_cpus" > input/ref_ILPs_temp_interpro.tsv || { echo "$(date) - InterProScan failed" >> pipeline.log; exit 1; }

# Count ILPs to balance with non-ILPs
ilp_count=$(grep -c "^>" input/ref_ILPs_temp.fasta)

# Fetch random non-ILP sequences to match ILP count
curl -o input/ref_non_ILPs.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)%20AND%20(taxonomy:%22Metazoa%20[33208]%22)%20AND%20NOT%20(insulin%20OR%20igf%20OR%20relaxin%20OR%20gonadulin%20OR%20bombyxin%20OR%20%22insulin-like%20peptide%22%20OR%20ilp%20OR%20neuroparsin)%20AND%20(length:[50%20TO%20200])&format=fasta&size=${ilp_count}&sort=rand" \
    || { echo "$(date) - UniProt non-ILP query failed" >> pipeline.log; exit 1; }

# Combine and annotate with Python script
cat input/ref_ILPs_temp.fasta input/ref_non_ILPs.fasta > input/ref_combined.fasta
python annotate_references.py input/ref_combined.fasta input/ref_ILPs_temp_interpro.tsv input/ref_ILPs.fasta || { echo "$(date) - annotate_references.py failed" >> pipeline.log; exit 1; }

# Clean up temporary files efficiently with pigz
pigz -f input/ref_ILPs_raw.fasta input/ref_ILPs_model.fasta input/ref_ILPs_nonmodel.fasta input/ref_non_ILPs.fasta input/ref_combined.fasta input/ref_ILPs_temp.fasta input/ref_ILPs_temp_interpro.tsv
echo "$(date) - 00a_fetch_references.sh completed" >> pipeline.log
touch input/.done
