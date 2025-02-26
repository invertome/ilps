#!/bin/bash
# Purpose: Step 04 - Generate structural models for identified ILP candidates
# Inputs: candidates/*_candidates.fasta
# Outputs: analysis/pdbs/*.pdb
# Config: config.yaml (max_cpus, colabfold_path)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Runs ColabFold on candidates only; includes checkpointing and input validation
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

# Process candidates into prepro, pro, mature forms and predict structures
for f in candidates/*_candidates.fasta; do
    base=$(basename "$f" _candidates.fasta)
    for type in prepro pro mature; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing $type for $base" >> pipeline.log
        for seq in $(seqkit seq -n "$f"); do
            fasta="analysis/pdbs/${type}_${base}_${seq}.fasta"
            pdb="analysis/pdbs/${type}_${base}_${seq}.pdb"
            python preprocess_ilp.py "$f" "$seq" "$type" "$fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type $seq" >> pipeline.log; exit 1; }
            if ! grep -q "$pdb" "$checkpoint" && [ ! -f "$pdb" ]; then
                "$colabfold_path" "$fasta" "analysis/pdbs/" --num-recycle 3 --model-type alphafold2_multimer_v3 --disulfide || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ColabFold failed for $fasta" >> pipeline.log; exit 1; }
                mv "analysis/pdbs/predict_*.pdb" "$pdb" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to move ColabFold output for $fasta" >> pipeline.log; exit 1; }
                echo "$pdb" >> "$checkpoint"
            fi
        done
    done
done

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 04_run_colabfold.sh completed in ${runtime}s" >> pipeline.log
touch analysis/.done_colabfold
