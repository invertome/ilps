#!/bin/bash
# Purpose: Preprocess input transcriptomes for ILP candidate identification
# Inputs: input/*.fasta (e.g., 9606_T1.fasta)
# Outputs: preprocess/*_preprocessed.fasta
# Config: config.yaml (max_cpus)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Translates nucleotide sequences, filters by length, deduplicates
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 01_preprocess.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f preprocess/.done_preprocess ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 01_preprocess.sh (already done)" >> pipeline.log
    exit 0
fi

# Check for input files
if ! ls input/[0-9]*_*.fasta >/dev/null 2>&1; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: No input FASTA files found in input/" >> pipeline.log
    exit 1
fi
if [ ! -f input/ref_ILPs.fasta ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: input/ref_ILPs.fasta not found" >> pipeline.log
    exit 1
fi

mkdir -p preprocess

python calc_ref_lengths.py input/ref_ILPs.fasta preprocess/ref_lengths.txt || { echo "$(date '+%Y-%m-%d %H:%M:%S') - calc_ref_lengths.py failed" >> pipeline.log; exit 1; }
min_len=$(awk 'NR==1{print $1-20}' preprocess/ref_lengths.txt)
max_len=$(awk 'NR==2{print $1+20}' preprocess/ref_lengths.txt)

for t in input/[0-9]*_*.fasta; do
    base=$(basename "$t" .fasta)
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing $t" >> pipeline.log
    seq_type=$(seqkit seq -t "$t" | head -n 1)
    if [ "$seq_type" == "DNA" ]; then
        TransDecoder.LongOrfs -t "$t" --threads "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - TransDecoder.LongOrfs failed for $t" >> pipeline.log; exit 1; }
        TransDecoder.Predict -t "$t" --cpu "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - TransDecoder.Predict failed for $t" >> pipeline.log; exit 1; }
        mv "${base}.transdecoder.pep" "preprocess/${base}_orfs.fasta"
        seqkit seq --threads "$max_cpus" -m "$min_len" -M "$max_len" "preprocess/${base}_orfs.fasta" > "preprocess/${base}_filtered.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - seqkit failed for $t" >> pipeline.log; exit 1; }
        cd-hit -i "preprocess/${base}_filtered.fasta" -o "preprocess/${base}_preprocessed.fasta" -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - CD-HIT failed for $t" >> pipeline.log; exit 1; }
        pigz -f "preprocess/${base}_orfs.fasta" "preprocess/${base}_filtered.fasta" "${base}.transdecoder"*
    else
        seqkit seq --threads "$max_cpus" -m "$min_len" -M "$max_len" "$t" > "preprocess/${base}_filtered.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - seqkit failed for $t" >> pipeline.log; exit 1; }
        cd-hit -i "preprocess/${base}_filtered.fasta" -o "preprocess/${base}_preprocessed.fasta" -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - CD-HIT failed for $t" >> pipeline.log; exit 1; }
        pigz -f "preprocess/${base}_filtered.fasta"
    fi
done

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 01_preprocess.sh completed in ${runtime}s" >> pipeline.log
touch preprocess/.done_preprocess
