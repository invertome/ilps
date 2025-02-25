#!/bin/bash
# Purpose: Identify ILP candidates from preprocessed transcriptomes using profile-based searches
# Inputs: preprocess/*_preprocessed.fasta (deduplicated protein sequences), input/ref_ILPs.fasta
# Outputs: 
#   - candidates/*_candidates.fasta (ILP candidates)
#   - Metadata: candidates/*_hhblits.out, *_hmm.out, *_hhsearch.out, *_blast.out, *_interpro.tsv
# Log: pipeline.log
# Notes: Uses linclust for clustering, batched BLAST for efficiency, parallel InterProScan
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f candidates/.done ]; then
    echo "$(date) - Skipping 02_identify_candidates.sh (already done)" >> pipeline.log
    exit 0
fi

# Set up directories
mkdir -p candidates hhblits_db
max_cpus=$(nproc)
echo "$(date) - Starting 02_identify_candidates.sh" >> pipeline.log

# Build HHblits database for candidate searches if not present
if [ ! -f hhblits_db/candidates_db ]; then
    echo "$(date) - Building HHblits database for candidates" >> pipeline.log
    cat preprocess/*_preprocessed.fasta > candidates/all_preprocessed.fasta
    hhblits -i candidates/all_preprocessed.fasta -d uniclav -o /dev/null -cpu "$max_cpus"
    mv uniclav* hhblits_db/
    ln -s hhblits_db/uniclav hhblits_db/candidates_db
fi

# Process each preprocessed file
for t in preprocess/[0-9]*_preprocessed.fasta; do
    base=$(basename "$t" _preprocessed.fasta)
    echo "Identifying candidates for $t..."
    # Cluster with linclust for redundancy reduction
    mmseqs linclust "$t" "candidates/${base}_clust" candidates/tmp --min-seq-id 0.95 --threads "$max_cpus" || { echo "$(date) - MMseqs2 linclust failed for $t" >> pipeline.log; exit 1; }
    # Search with HHblits and HMMER for ILP similarity
    hhblits -i "candidates/${base}_clust_rep_seq.fasta" -d hhblits_db/candidates_db -o "candidates/${base}_hhblits.out" -n 3 -cpu "$max_cpus" -v 1 || { echo "$(date) - HHblits failed for $t" >> pipeline.log; exit 1; }
    hmmsearch -o "candidates/${base}_hmm.out" --cpu "$max_cpus" input/ilp.hmm "candidates/${base}_clust_rep_seq.fasta" || { echo "$(date) - HMMER failed for $t" >> pipeline.log; exit 1; }
    python combine_hits.py "candidates/${base}_hhblits.out" "candidates/${base}_hmm.out" "candidates/${base}_candidates.fasta" || { echo "$(date) - combine_hits.py failed for $t" >> pipeline.log; exit 1; }
    # Build HMM profile and search
    hhmake -i "candidates/${base}_candidates.fasta" -o "candidates/${base}_candidates.hhm" -v 0 || { echo "$(date) - hhmake failed for $t" >> pipeline.log; exit 1; }
    hhsearch -i "candidates/${base}_candidates.hhm" -d input/ilp_db.hhm -o "candidates/${base}_hhsearch.out" -cpu "$max_cpus" -v 1 || { echo "$(date) - HHsearch failed for $t" >> pipeline.log; exit 1; }
    pigz -f candidates/${base}_clust*
done

# Batch BLAST against reference ILPs for identity scoring
makeblastdb -in input/ref_ILPs.fasta -dbtype prot -out input/ref_ILPs_db || { echo "$(date) - makeblastdb failed" >> pipeline.log; exit 1; }
blastp -query <(cat candidates/*_candidates.fasta) -db input/ref_ILPs_db -out candidates/all_blast.out -outfmt 6 -num_threads "$max_cpus" -evalue 1e-5 || { echo "$(date) - BLAST failed" >> pipeline.log; exit 1; }
for f in candidates/*_candidates.fasta; do
    base=$(basename "$f" _candidates.fasta)
    grep "^${base}_" candidates/all_blast.out > "candidates/${base}_blast.out"
done

# Parallel InterProScan for domain annotation
ls candidates/*_candidates.fasta | parallel -j "$max_cpus" "interproscan.sh -i {} -dp -f tsv -iprlookup -goterms -pa -cpu 1 > candidates/{/.}_interpro.tsv" || { echo "$(date) - InterProScan failed" >> pipeline.log; exit 1; }

# Clean up temporary files
pigz -f candidates/tmp/*
echo "$(date) - 02_identify_candidates.sh completed" >> pipeline.log
touch candidates/.done
