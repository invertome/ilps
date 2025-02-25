#!/bin/bash
# Purpose: Run ColabFold to generate structural models for all sequence types
# Inputs: preprocess/*_${type}.fasta, analysis/pdbs/*_${type}.fasta
# Outputs: preprocess/*.pdb, analysis/pdbs/*.pdb
# Config: config.yaml (max_cpus, colabfold_path)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Includes checkpointing for error recovery
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
colabfold_path=$(yq e '.colabfold_path' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 00c_run_colabfold.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f colabfold.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 00c_run_colabfold.sh (already done)" >> pipeline.log
    exit 0
fi

# Checkpoint file
checkpoint="colabfold_checkpoint.txt"
touch "$checkpoint"

for type in prepro pro mature; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Modeling $type sequences" >> pipeline.log
    # Process reference sequences
    ls preprocess/*_${type}.fasta | while read -r fasta; do
        pdb="preprocess/$(basename "$fasta" .fasta).pdb"
        if ! grep -q "$pdb" "$checkpoint" && [ ! -f "$pdb" ]; then
            "$colabfold_path" "$fasta" "$(dirname "$pdb")" --num-recycle 3 --model-type alphafold2_multimer_v3 --disulfide || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ColabFold failed for $fasta" >> pipeline.log; exit 1; }
            mv "$(dirname "$pdb")/predict_*.pdb" "$pdb"
            echo "$pdb" >> "$checkpoint"
        fi
    done
    # Process candidate sequences if directory exists
    if [ -d analysis/pdbs ]; then
        ls analysis/pdbs/*_${type}.fasta | while read -r fasta; do
            pdb="analysis/pdbs/$(basename "$fasta" .fasta).pdb"
            if ! grep -q "$pdb" "$checkpoint" && [ ! -f "$pdb" ]; then
                "$colabfold_path" "$fasta" "$(dirname "$pdb")" --num-recycle 3 --model-type alphafold2_multimer_v3 --disulfide || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ColabFold failed for $fasta" >> pipeline.log; exit 1; }
                mv "$(dirname "$pdb")/predict_*.pdb" "$pdb"
                echo "$pdb" >> "$checkpoint"
            fi
        done
    fi
done

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 00c_run_colabfold.sh completed in ${runtime}s" >> pipeline.log
touch colabfold.done
