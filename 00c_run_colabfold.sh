#!/bin/bash
# Purpose: Run ColabFold to generate structural models for prepropeptides, propeptides, and mature peptides
# Inputs: preprocess/*_${type}.fasta, analysis/pdbs/*_${type}.fasta (type: prepro, pro, mature)
# Outputs: preprocess/*.pdb, analysis/pdbs/*.pdb (structural models for all types)
# Log: pipeline.log
# Notes: Optimized for GPU execution; skips existing PDBs for incremental processing
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f colabfold.done ]; then
    echo "$(date) - Skipping 00c_run_colabfold.sh (already done)" >> pipeline.log
    exit 0
fi

max_cpus=$(nproc)
echo "$(date) - Starting 00c_run_colabfold.sh" >> pipeline.log

# Process each sequence type (prepro, pro, mature) in parallel, skipping existing PDBs
for type in prepro pro mature; do
    echo "$(date) - Modeling $type sequences" >> pipeline.log
    # Generate models for reference sequences if PDBs don't exist
    ls preprocess/*_${type}.fasta | parallel -j "$max_cpus" "[ -f preprocess/{/.}.pdb ] || python run_alphafold2.py {} preprocess/{/.}.pdb --disulfide_constraints --num_recycle 3" || { echo "$(date) - ColabFold failed for preprocess $type" >> pipeline.log; exit 1; }
    # Generate models for candidate sequences if analysis directory exists
    if [ -d analysis/pdbs ]; then
        ls analysis/pdbs/*_${type}.fasta | parallel -j "$max_cpus" "[ -f analysis/pdbs/{/.}.pdb ] || python run_alphafold2.py {} analysis/pdbs/{/.}.pdb --disulfide_constraints --num_recycle 3" || { echo "$(date) - ColabFold failed for analysis $type" >> pipeline.log; exit 1; }
    fi
done

echo "$(date) - 00c_run_colabfold.sh completed" >> pipeline.log
touch colabfold.done
