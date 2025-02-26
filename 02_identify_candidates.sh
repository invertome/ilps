#!/bin/bash
# Purpose: Step 02 - Identify ILP candidates from preprocessed transcriptomes
# Inputs: preprocess/*_preprocessed.fasta, input/ref_ILPs.fasta, input/ilp.hmm, input/ilp_db.hhm
# Outputs: candidates/*_candidates.fasta, metadata files
# Config: config.yaml (max_cpus, interpro_path)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Caches BLAST and InterProScan results; added input checks
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

# Check required input files
for file in preprocess/[0-9]*_preprocessed.fasta input/ref_ILPs.fasta input/ilp.hmm input/ilp_db.hhm; do
    if ! ls "$file" >/dev/null 2>&1; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: $file not found" >> pipeline.log
        exit 1
    fi
done

mkdir -p candidates hhblits_db

if [ ! -f hhblits_db/candidates_db ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Building HHblits database for candidates" >> pipeline.log
    cat preprocess/*_preprocessed.fasta > candidates/all_preprocessed.fasta
    hhblits -i candidates/all_preprocessed.fasta -d uniclav -o /dev/null -cpu "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHblits database build failed" >> pipeline.log; exit 1; }
    mv uniclav* hhblits_db/
    ln -s hhblits_db/uniclav hhblits_db/candidates_db
fi

for t in preprocess/[0-9]*_preprocessed.fasta; do
    base=$(basename "$t" _preprocessed.fasta)
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Identifying candidates for $t" >> pipeline.log
    mmseqs linclust "$t" "candidates/${base}_clust" candidates/tmp --min-seq-id 0.95 --threads "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 linclust failed for $t" >> pipeline.log; exit 1; }
    hhblits -i "candidates/${base}_clust_rep_seq.fasta" -d hhblits_db/candidates_db -o "candidates/${base}_hhblits.out" -n 3 -cpu "$max_cpus" -v 1 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHblits failed for $t" >> pipeline.log; exit 1; }
    hmmsearch -o "candidates/${base}_hmm.out" --cpu "$max_cpus" input/ilp.hmm "candidates/${base}_clust_rep_seq.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HMMER failed for $t" >> pipeline.log; exit 1; }
    python combine_hits.py "candidates/${base}_hhblits.out" "candidates/${base}_hmm.out" "candidates/${base}_candidates.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - combine_hits.py failed for $t" >> pipeline.log; exit 1; }
    hhmake -i "candidates/${base}_candidates.fasta" -o "candidates/${base}_candidates.hhm" -v 0 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - hhmake failed for $t" >> pipeline.log; exit 1; }
    hhsearch -i "candidates/${base}_candidates.hhm" -d input/ilp_db.hhm -o "candidates/${base}_hhsearch.out" -cpu "$max_cpus" -v 1 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHsearch failed for $t" >> pipeline.log; exit 1; }
    pigz -f candidates/${base}_clust*
done

if [ ! -f candidates/all_blast.out ]; then
    makeblastdb -in input/ref_ILPs.fasta -dbtype prot -out input/ref_ILPs_db || { echo "$(date '+%Y-%m-%d %H:%M:%S') - makeblastdb failed" >> pipeline.log; exit 1; }
    blastp -query <(cat candidates/*_candidates.fasta) -db input/ref_ILPs_db -out candidates/all_blast.out -outfmt 6 -num_threads "$max_cpus" -evalue 1e-5 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - BLAST failed" >> pipeline.log; exit 1; }
fi
for f in candidates/*_candidates.fasta; do
    base=$(basename "$f" _candidates.fasta)
    grep "^${base}_" candidates/all_blast.out > "candidates/${base}_blast.out"
done

ls candidates/*_candidates.fasta | while read -r f; do
    base=$(basename "$f" .fasta)
    if [ ! -f "candidates/${base}_interpro.tsv" ]; then
        "$interpro_path" -i "$f" -dp -f tsv -iprlookup -goterms -pa -cpu 1 > "candidates/${base}_interpro.tsv" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - InterProScan failed for $f" >> pipeline.log; exit 1; }
    fi
done

pigz -f candidates/tmp/*
end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 02_identify_candidates.sh completed in ${runtime}s" >> pipeline.log
touch candidates/.done
