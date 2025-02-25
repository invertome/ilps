#!/bin/bash
# Purpose: Step 05 - Perform comparative and phylogenetic analysis on ILP candidates
# Inputs: analysis/all_candidates.fasta, analysis/predictions.csv, analysis/novel_candidates.csv
# Outputs: analysis/ (trees, pdbs/, clades_*, unaligned_*, aligned_*, plots_*)
# Config: config.yaml (max_cpus)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Generates all possible consensus trees from sequence and structural data; ensures environment consistency
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
cpus_per_task=$((max_cpus / 3))
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 05_comparative_analysis.sh with $cpus_per_task CPUs per task" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f analysis/.done_comparative ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 05_comparative_analysis.sh (already done)" >> pipeline.log
    exit 0
fi

# Check dependencies
if ! command -v mamba >/dev/null 2>&1 || ! command -v snakemake >/dev/null 2>&1; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: mamba or snakemake not installed" >> pipeline.log
    exit 1
fi

mkdir -p analysis/pdbs clades_ete_prepro clades_autophy_prepro clades_ete_pro clades_autophy_pro clades_ete_mature clades_autophy_mature analysis/unaligned_prepro analysis/aligned_prepro analysis/plots_prepro analysis/unaligned_pro analysis/aligned_pro analysis/plots_pro analysis/unaligned_mature analysis/aligned_mature analysis/plots_mature

python filter_ilps.py analysis/all_candidates.fasta analysis/predictions.csv analysis/ilp_candidates.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - filter_ilps.py failed" >> pipeline.log; exit 1; }
python determine_common_taxonomy.py input/[0-9]*_*.fasta analysis/common_taxonomy.txt || { echo "$(date '+%Y-%m-%d %H:%M:%S') - determine_common_taxonomy.py failed" >> pipeline.log; exit 1; }
common_rank=$(cat analysis/common_taxonomy.txt | cut -d' ' -f1)
common_taxid=$(cat analysis/common_taxonomy.txt | cut -d' ' -f2)
python filter_ref_ilps_by_taxonomy.py input/ref_ILPs.fasta "$common_taxid" analysis/ref_ILPs_filtered.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - filter_ref_ilps_by_taxonomy.py failed" >> pipeline.log; exit 1; }

for type in prepro pro mature; do
    for seq in $(seqkit seq -n analysis/ref_ILPs_filtered.fasta); do
        python preprocess_ilp.py analysis/ref_ILPs_filtered.fasta "$seq" "$type" "preprocess/${type}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type (ref)" >> pipeline.log; exit 1; }
    done
    for seq in $(seqkit seq -n analysis/ilp_candidates.fasta); do
        python preprocess_ilp.py analysis/ilp_candidates.fasta "$seq" "$type" "analysis/pdbs/${type}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type (candidate)" >> pipeline.log; exit 1; }
    done
done

echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing sequence alignments and trees" >> pipeline.log
checkpoint="iqtree_checkpoint.txt"
touch "$checkpoint"
echo "prepro pro mature" | parallel -j 3 --eta "\
    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting processing for {}\"; \
    cat analysis/ilp_candidates.fasta analysis/ref_ILPs_filtered.fasta > analysis/{}_analysis_candidates.fasta || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - cat failed for {}\"; exit 1; }; \
    if [ ! -f analysis/{}_alignment_trimmed.fasta ]; then \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running MAFFT for {}\"; \
        mafft --globalpair --maxiterate 1000 --thread $cpus_per_task analysis/{}_analysis_candidates.fasta > analysis/{}_alignment.fasta || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - MAFFT failed for {}\"; exit 1; }; \
        trimal -in analysis/{}_alignment.fasta -out analysis/{}_alignment_trimmed.fasta -automated1 || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - trimal failed for {}\"; exit 1; }; \
    fi; \
    FastTree -lg -nt -fastest analysis/{}_alignment_trimmed.fasta > analysis/{}_fasttree.tre || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - FastTree failed for {}\"; exit 1; }; \
    if ! grep -q \"analysis/{}_iqtree.treefile\" \"$checkpoint\"; then \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running IQ-TREE for {}\"; \
        iqtree -s analysis/{}_alignment_trimmed.fasta -m TEST -abayes -B 10000 -bnni -T $cpus_per_task -t analysis/{}_fasttree.tre -pre analysis/{}_iqtree --redo || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - IQ-TREE failed for {}\"; exit 1; }; \
        echo \"analysis/{}_iqtree.treefile\" >> \"$checkpoint\"; \
    fi; \
    mv analysis/{}_iqtree.treefile analysis/{}_alignment.fasta.treefile"

# Generate structural trees using Foldtree Snakemake workflow
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating structural trees with Foldtree" >> pipeline.log

# Ensure running in ilp_pipeline environment
if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "ilp_pipeline" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Warning: Not in ilp_pipeline environment; activating it" >> pipeline.log
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ilp_pipeline || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to activate ilp_pipeline" >> pipeline.log; exit 1; }
fi

if [ ! -d "fold_tree" ]; then
