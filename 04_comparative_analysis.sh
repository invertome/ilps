#!/bin/bash
# Purpose: Perform comparative and phylogenetic analysis on ILP candidates for prepropeptides, propeptides, and mature peptides
# Inputs: 
#   - analysis/all_candidates.fasta (all candidates)
#   - analysis/predictions.csv (ML probabilities)
#   - analysis/novel_candidates.csv (novelty flags)
# Outputs: 
#   - analysis/ (trees, pdbs/, clades_ete_*/, clades_autophy_*/, unaligned_*/aligned_*/plots_*)
#   - Multiple consensus trees (sequence and structural combinations for all types)
# Log: pipeline.log
# Notes: Builds models and structural trees for all types; includes trimming, parallel structural tree generation, and optimized CPU allocation
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f analysis/.done_comparative ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 04_comparative_analysis.sh (already done)" >> pipeline.log
    exit 0
fi

# Set up directories for each sequence type
mkdir -p analysis/pdbs clades_ete_prepro clades_autophy_prepro clades_ete_pro clades_autophy_pro clades_ete_mature clades_autophy_mature analysis/unaligned_prepro analysis/aligned_prepro analysis/plots_prepro analysis/unaligned_pro analysis/aligned_pro analysis/plots_pro analysis/unaligned_mature analysis/aligned_mature analysis/plots_mature
max_cpus=$(nproc)
cpus_per_task=$((max_cpus / 3))  # Allocate CPUs per task, limiting to 3 concurrent runs for efficiency
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 04_comparative_analysis.sh with $cpus_per_task CPUs per task" >> pipeline.log

# Filter candidates to ILPs based on ML probability (> 0.7)
python filter_ilps.py analysis/all_candidates.fasta analysis/predictions.csv analysis/ilp_candidates.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - filter_ilps.py failed" >> pipeline.log; exit 1; }

# Determine highest common taxonomic rank for reference filtering
python determine_common_taxonomy.py input/[0-9]*_*.fasta analysis/common_taxonomy.txt || { echo "$(date '+%Y-%m-%d %H:%M:%S') - determine_common_taxonomy.py failed" >> pipeline.log; exit 1; }
common_rank=$(cat analysis/common_taxonomy.txt | cut -d' ' -f1)
common_taxid=$(cat analysis/common_taxonomy.txt | cut -d' ' -f2)

# Filter reference ILPs to common taxonomy and prepare for modeling
python filter_ref_ilps_by_taxonomy.py input/ref_ILPs.fasta "$common_taxid" analysis/ref_ILPs_filtered.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - filter_ref_ilps_by_taxonomy.py failed" >> pipeline.log; exit 1; }

# Generate FASTA files for all sequence types (prepro, pro, mature) for references and candidates
mkdir -p analysis/pdbs
for type in prepro pro mature; do
    for seq in $(seqkit seq -n analysis/ref_ILPs_filtered.fasta); do
        python preprocess_ilp.py analysis/ref_ILPs_filtered.fasta "$seq" "$type" "preprocess/${type}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type (reference)" >> pipeline.log; exit 1; }
    done
    for seq in $(seqkit seq -n analysis/ilp_candidates.fasta); do
        python preprocess_ilp.py analysis/ilp_candidates.fasta "$seq" "$type" "analysis/pdbs/${type}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed for $type (candidate)" >> pipeline.log; exit 1; }
    done
done

# Process each sequence type (prepro, pro, mature) with optimized parallelization
echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing sequence alignments and trees" >> pipeline.log
echo "prepro pro mature" | parallel -j 3 --eta "\
    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Starting processing for {}\"; \
    cat analysis/ilp_candidates.fasta analysis/ref_ILPs_filtered.fasta > analysis/{}_analysis_candidates.fasta || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - cat failed for {}\"; exit 1; }; \
    if [ ! -f analysis/{}_alignment_trimmed.fasta ]; then \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running MAFFT for {}\"; \
        mafft --globalpair --maxiterate 1000 --thread $cpus_per_task analysis/{}_analysis_candidates.fasta > analysis/{}_alignment.fasta || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - MAFFT failed for {}\"; exit 1; }; \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running trimal for {}\"; \
        trimal -in analysis/{}_alignment.fasta -out analysis/{}_alignment_trimmed.fasta -automated1 || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - trimal failed for {}\"; exit 1; }; \
    fi; \
    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running FastTree for {}\"; \
    FastTree -lg -nt -fastest analysis/{}_alignment_trimmed.fasta > analysis/{}_fasttree.tre || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - FastTree failed for {}\"; exit 1; }; \
    echo \"$(date '+%Y-%m-%d %H:%M:%S') - Running IQ-TREE for {}\"; \
    iqtree -s analysis/{}_alignment_trimmed.fasta -m TEST -abayes -B 10000 -bnni -T $cpus_per_task -t analysis/{}_fasttree.tre -pre analysis/{}_iqtree --redo || { echo \"$(date '+%Y-%m-%d %H:%M:%S') - IQ-TREE failed for {}\"; exit 1; }; \
    mv analysis/{}_iqtree.treefile analysis/{}_alignment.fasta.treefile"

# Generate structural trees in parallel for all types with fewer concurrent runs
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating structural trees" >> pipeline.log
echo "prepro pro mature" | parallel -j 3 --eta "foldtree -i analysis/pdbs/ -o analysis/{}_structural_tree.tre -type {}" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - FoldTree failed" >> pipeline.log; exit 1; }

# Generate consensus trees for each type
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating per-type consensus trees" >> pipeline.log
for type in prepro pro mature; do
    Rscript consensus_tree.R "analysis/${type}_alignment.fasta.treefile" "analysis/${type}_structural_tree.tre" "analysis/${type}_consensus_tree.tre" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for $type" >> pipeline.log; exit 1; }
done

# Generate additional consensus trees for sequence combinations
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating additional sequence consensus trees" >> pipeline.log
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/pro_alignment.fasta.treefile analysis/prepro_pro_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + pro (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/pro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/pro_mature_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for pro + mature (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/prepro_mature_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + mature (sequence)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_alignment.fasta.treefile analysis/pro_alignment.fasta.treefile analysis/mature_alignment.fasta.treefile analysis/all_seq_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for all three (sequence)" >> pipeline.log; exit 1; }

# Generate additional consensus trees for structural combinations
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating additional structural consensus trees" >> pipeline.log
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/pro_structural_tree.tre analysis/prepro_pro_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + pro (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/pro_structural_tree.tre analysis/mature_structural_tree.tre analysis/pro_mature_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for pro + mature (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/mature_structural_tree.tre analysis/prepro_mature_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for prepro + mature (structural)" >> pipeline.log; exit 1; }
Rscript consensus_tree.R analysis/prepro_structural_tree.tre analysis/pro_structural_tree.tre analysis/mature_structural_tree.tre analysis/all_struct_consensus.tre || { echo "$(date '+%Y-%m-%d %H:%M:%S') - consensus_tree.R failed for all three (structural)" >> pipeline.log; exit 1; }

# Clade analysis with ETE and Autophy for each type
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

echo "$(date '+%Y-%m-%d %H:%M:%S') - 04_comparative_analysis.sh completed" >> pipeline.log
touch analysis/.done_comparative
