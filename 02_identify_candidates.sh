#!/bin/bash
# Purpose: Identify ILP candidates from preprocessed transcriptomes
# Inputs: preprocess/*_preprocessed.fasta, input/ref_ILPs.fasta
# Outputs: candidates/*_candidates.fasta, metadata files
# Config: config.yaml (max_cpus, interpro_path)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Caches BLAST and InterProScan results
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
interpro_path=$(yq e '.interpro_path' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 02_identify_candidates.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f candidates/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 02_identify_candidates.sh (already done)" >> pipeline.log
    exit 0
fi

mkdir -p candidates hhblits_db

# Build HHblits database
if [ ! -f hhblits_db/candidates_db ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Building HHblits database for candidates" >> pipeline.log
    cat preprocess/*_preprocessed.fasta > candidates/all_preprocessed.fasta
    hhblits -i candidates/all_preprocessed.fasta -d uniclav -o /dev/null -cpu "$max_cpus"
    mv uniclav* hhblits_db/
    ln -s hhblits_db/uniclav hhblits_db/candidates_db
fi

for t in preprocess/[0-9]*_preprocessed.fasta; do
    base=$(basename "$t" _preprocessed.fasta)
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Identifying candidates for $t" >> pipeline.log
    mmseqs linclust "$t" "candidates/${base}_clust" candidates/tmp --min-seq-id 0.95 --threads "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 linclust failed for $t" >> pipeline.log; exit 1; }
    hhblits -i "candidates/${base}_clust_rep_seq.fasta" -d hhblits_db/candidates_db -o "candidates/${base}_hhblits.out" -n 3 -cpu "$max_cpus" -v 1 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHblits failed for $t" >> pipeline.log; exit 1; }
    hmmsearch -o "candidates/${base}_hmm.out" --cpu "$max_cpus" input/ilp.hmm "candidates/${base}_clust_rep_seq.fasta
