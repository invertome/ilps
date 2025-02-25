#!/bin/bash
# Purpose: Preprocess input transcriptomes (nucleotide or amino acid) for ILP candidate identification
# Inputs: input/*.fasta (e.g., 9606_T1.fasta)
# Outputs: preprocess/*_preprocessed.fasta (deduplicated protein sequences)
# Log: pipeline.log
# Notes: Translates nucleotide sequences with TransDecoder, filters by length based on reference ILPs, deduplicates at 95%
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f preprocess/.done_preprocess ]; then
    echo "$(date) - Skipping 01_preprocess.sh (already done)" >> pipeline.log
    exit 0
fi

mkdir -p preprocess
max_cpus=$(nproc)
echo "$(date) - Starting 01_preprocess.sh" >> pipeline.log

# Calculate dynamic length thresholds from reference ILPs
python calc_ref_lengths.py input/ref_ILPs.fasta preprocess/ref_lengths.txt
if [ $? -ne 0 ]; then echo "$(date) - Error in calc_ref_lengths.py" >> pipeline.log; exit 1; fi
min_len=$(awk 'NR==1{print $1-20}' preprocess/ref_lengths.txt)  # Min length minus buffer
max_len=$(awk 'NR==2{print $1+20}' preprocess/ref_lengths.txt)  # Max length plus buffer

# Process each input FASTA file
for t in input/[0-9]*_*.fasta; do
    base=$(basename "$t" .fasta)
    echo "Processing $t..."
    seq_type=$(seqkit seq -t "$t" | head -n 1)
    if [ "$seq_type" == "DNA" ]; then
        # Translate nucleotide sequences and extract ORFs
        TransDecoder.LongOrfs -t "$t" --threads "$max_cpus" || { echo "$(date) - TransDecoder.LongOrfs failed for $t" >> pipeline.log; exit 1; }
        TransDecoder.Predict -t "$t" --cpu "$max_cpus" || { echo "$(date) - TransDecoder.Predict failed for $t" >> pipeline.log; exit 1; }
        mv "${base}.transdecoder.pep" "preprocess/${base}_orfs.fasta"
        # Filter by length and deduplicate
        seqkit seq --threads "$max_cpus" -m "$min_len" -M "$max_len" "preprocess/${base}_orfs.fasta" > "preprocess/${base}_filtered.fasta" || { echo "$(date) - seqkit failed for $t" >> pipeline.log; exit 1; }
        cd-hit -i "preprocess/${base}_filtered.fasta" -o "preprocess/${base}_preprocessed.fasta" -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date) - CD-HIT failed for $t" >> pipeline.log; exit 1; }
        pigz -f "preprocess/${base}_orfs.fasta" "preprocess/${base}_filtered.fasta" "${base}.transdecoder"*
    else
        # Filter amino acid sequences by length and deduplicate
        seqkit seq --threads "$max_cpus" -m "$min_len" -M "$max_len" "$t" > "preprocess/${base}_filtered.fasta" || { echo "$(date) - seqkit failed for $t" >> pipeline.log; exit 1; }
        cd-hit -i "preprocess/${base}_filtered.fasta" -o "preprocess/${base}_preprocessed.fasta" -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date) - CD-HIT failed for $t" >> pipeline.log; exit 1; }
        pigz -f "preprocess/${base}_filtered.fasta"
    fi
done

echo "$(date) - 01_preprocess.sh completed" >> pipeline.log
touch preprocess/.done_preprocess
