#!/bin/bash
# Purpose: Step 04 - Generate structural models for identified ILP candidates using ColabFold
# Inputs: candidates/*_candidates.fasta
# Outputs: analysis/pdbs/*.pdb
# Config: config.yaml (max_cpus, colabfold_path)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Uses DeepNeuropePred via preprocess_ilp.py; parallelized with GNU Parallel
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
colabfold_path=$(yq e '.colabfold_path' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 04_run_colabfold.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f analysis/.done_colabfold ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 04_run_colabfold.sh (already done)" >> pipeline.log
    exit 0
fi

# Check for input candidate FASTA files
if ! ls candidates/*_candidates.fasta >/dev/null 2>&1; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: No candidate FASTA files found in candidates/" >> pipeline.log
    exit 1
fi

mkdir -p analysis/pdbs
checkpoint="analysis/colabfold_checkpoint.txt"
touch "$checkpoint"

# Function to process a single sequence
process_sequence() {
    fasta="$1"
    pdb="$2"
    type="$3"
    base="$4"
    seq="$5"

    if ! grep -q "$pdb" "$checkpoint" && [ ! -f "$pdb" ]; then
        python preprocess_ilp.py "$fasta" "$seq" "$type" "analysis/pdbs/${type}_${base}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type $seq" >> pipeline.log; exit 1; }
        "$colabfold_path" "analysis/pdbs/${type}_${base}_${seq}.fasta" "analysis/pdbs/" --num-recycle 3 --model-type alphafold2_multimer_v3 --disulfide || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ColabFold failed for $seq" >> pipeline.log; exit 1; }
        mv "analysis/pdbs/predict_*.pdb" "$pdb" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to move ColabFold output for $seq" >> pipeline.log; exit 1; }
        echo "$pdb" >> "$checkpoint"
    fi
}

export -f process_sequence
export checkpoint
export colabfold_path

# Generate jobs for parallel processing
rm -f colabfold_jobs.txt
for f in candidates/*_candidates.fasta; do
    base=$(basename "$f" _candidates.fasta)
    for type in prepro pro mature; do
        seq_ids=$(seqkit seq -n "$f")
        for seq in $seq_ids; do
            echo "process_sequence $f analysis/pdbs/${type}_${base}_${seq}.pdb $type $base $seq" >> colabfold_jobs.txt
        done
    done
done

parallel -j "$max_cpus" < colabfold_jobs.txt || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Parallel ColabFold failed" >> pipeline.log; exit 1; }

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 04_run_colabfold.sh completed in ${runtime}s" >> pipeline.log
touch analysis/.done_colabfold
rm -f colabfold_jobs.txt
