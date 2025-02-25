#!/bin/bash
# Purpose: Perform comparative and phylogenetic analysis on ILP candidates
# Inputs: analysis/all_candidates.fasta, analysis/predictions.csv, analysis/novel_candidates.csv
# Outputs: analysis/ (trees, pdbs/, clades_*, unaligned_*, aligned_*, plots_*)
# Config: config.yaml (max_cpus)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Checkpointing for IQ-TREE, optimized CPU allocation
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
cpus_per_task=$((max_cpus / 3))
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 04_comparative_analysis.sh with $cpus_per_task CPUs per task" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f analysis/.done_comparative ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 04_comparative_analysis.sh (already done)" >> pipeline.log
    exit 0
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

# Ensure Foldtree repository is cloned and environment is set up
if [ ! -d "fold_tree" ]; then
    git clone git@github.com:DessimozLab/fold_tree.git || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to clone Foldtree repository" >> pipeline.log; exit 1; }
    mamba env create -n foldtree --file=./fold_tree/workflow/config/fold_tree.yaml || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to create Foldtree environment" >> pipeline.log; exit 1; }
fi

# Activate Foldtree environment
source "$(conda info --base)/etc/profile.d/conda.sh"
mamba activate foldtree || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to activate Foldtree environment" >> pipeline.log; exit 1; }

# Process each sequence type
for type in prepro pro mature; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Setting up Foldtree for $type" >> pipeline.log
    
    # Create working directory and structs subdirectory
    foldtree_dir="analysis/${type}_foldtree"
    mkdir -p "$foldtree_dir/structs" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to create $foldtree_dir/structs" >> pipeline.log; exit 1; }
    
    # Copy PDB files to structs directory
    cp analysis/pdbs/${type}_*.pdb "$foldtree_dir/structs/" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to copy PDBs for $type" >> pipeline.log; exit 1; }
    
    # Create an empty identifiers file (required even for custom structs)
    touch "$foldtree_dir/identifiers.txt"
    
    # Run Foldtree Snakemake workflow
    snakemake --cores "$cpus_per_task" --use-conda -s ./fold_tree/workflow/fold_tree \
        --config folder="$foldtree_dir" filter=False custom_structs=True foldseek_cores="$cpus_per_task" \
        || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Foldtree failed for $type" >> pipeline.log; exit 1; }
    
    # Move the Foldtree output tree to the expected location
    mv "$foldtree_dir/Foldtree.tre" "analysis/${type}_structural_tree.tre" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to move Foldtree output for $type" >> pipeline.log; exit 1; }
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Foldtree completed for $type" >> pipeline.log
done

# Deactivate Foldtree environment
mamba deactivate

for type in prepro pro mature; do
    Rscript consensus_tree.R "analysis/${type}_alignment.fasta.treefile" "analysis/${type}_structural_tree.tre" "analysis/${type}_consensus_tree.tre" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for $type" >> pipeline.log; exit 1; }
done

# Additional consensus trees (sequence)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating additional sequence consensus trees" >> pipeline.log
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/pro_alignment.fasta.treefile analysis/prepro_pro_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + pro (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/pro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/pro_mature_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for pro + mature (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/prepro_mature_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + mature (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/pro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/all_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for all three (sequence)" >> pipeline.log; exit 1; }

# Additional consensus trees (structural)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating additional structural consensus trees" >> pipeline.log
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/pro_structural_tree.tre analysis/prepro_pro_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + pro (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/pro_structural_tree.tre analysis/mature_structural_tree.tre analysis/pro_mature_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for pro + mature (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/mature_structural_tree.tre analysis/prepro_mature_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + mature (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/pro_structural_tree.tre analysis/mature_structural_tree.tre analysis/all_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for all three (structural)" >> pipeline.log; exit 1; }

for type in prepro pro mature; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting clade analysis for $type" >> pipeline.log
    python identify_clades.py "analysis/${type}_alignment.fasta.treefile" "clades_ete_${type}/" "analysis/unaligned_${type}/" "analysis/aligned_${type}/" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - identify_clades.py failed for $type" >> pipeline.log; exit 1; }
    autophy -t "analysis/${type}_alignment.fasta.treefile" -o "clades_autophy_${type}/" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Autophy failed for $type" >> pipeline.log; exit 1; }
    for clade_dir in "clades_ete_${type}" "clades_autophy_${type}"; do
        method=$(basename "$clade_dir" | cut -d'_' -f2)
        ls "$clade_dir"/*.fasta | parallel -j "$max_cpus" "meme {} -o ${clade_dir}/meme_{/} -maxw 20 -nmotifs 5 -mod zoops" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - MEME failed for $method $type" >> pipeline.log; exit 1; }
        ls "$clade_dir"/*.fasta | parallel -j "$max_cpus" "ame --control analysis/${type}_analysis_candidates.fasta ${clade_dir}/meme_{/}/meme.txt > ${clade_dir}/ame_{/}.txt" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - AME failed for $method $type" >> pipeline.log; exit 1; }
        ls "$clade_dir"/*.fasta | parallel -j "$max_cpus" "fimo ${clade_dir}/meme_{/}/meme.txt analysis/${type}_analysis_candidates.fasta > ${clade_dir}/fimo_{/}.txt" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - FIMO failed for $method $type" >> pipeline.log; exit 1; }
        ls "$clade_dir"/*.fasta | parallel -j "$max_cpus" "mafft --globalpair --maxiterate 1000 --thread 1 {} > analysis/aligned_${type}/${method}_{/}_aligned.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - MAFFT clade alignment failed for $method $type" >> pipeline.log; exit 1; }
        ls "analysis/aligned_${type}/${method}"/*_aligned.fasta | parallel -j "$max_cpus" "python plot_alignment.py {} analysis/plots_${type}/${method}_{/}_logo.png" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - plot_alignment.py failed for $method $type" >> pipeline.log; exit 1; }
    done
done

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 04_comparative_analysis.sh completed in ${runtime}s" >> pipeline.log
touch analysis/.done_comparative
