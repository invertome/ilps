#!/bin/bash
# Purpose: Prepare training data for machine learning by processing reference sequences
# Inputs: input/ref_ILPs.fasta (annotated ILP/non-ILP sequences)
# Outputs: 
#   - preprocess/ref_features.csv (feature matrix for ML)
#   - preprocess/ref_labels.csv (labels: 1 for ILP, 0 for non-ILP)
#   - preprocess/*_${type}.fasta (prepro, pro, mature sequences for modeling)
# Log: pipeline.log
# Notes: Uses linclust for clustering, SignalP for signal peptide removal, and InterPro for domain annotation
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

if [ -f preprocess/.done ]; then
    echo "$(date) - Skipping 00b_prepare_training.sh (already done)" >> pipeline.log
    exit 0
fi

# Set up directories
mkdir -p preprocess hhblits_db
max_cpus=$(nproc)
echo "$(date) - Starting 00b_prepare_training.sh" >> pipeline.log

# Build HHblits database if not already present
if [ ! -f hhblits_db/ref_db ]; then
    echo "$(date) - Building HHblits database from UniClust" >> pipeline.log
    hhblits -i input/ref_ILPs.fasta -d uniclav -o /dev/null -cpu "$max_cpus"
    mv uniclav* hhblits_db/
    ln -s hhblits_db/uniclav hhblits_db/ref_db
fi

# Cluster sequences with linclust for redundancy reduction
mmseqs linclust input/ref_ILPs.fasta preprocess/ref_clust preprocess/tmp --min-seq-id 0.95 --threads "$max_cpus" || { echo "$(date) - MMseqs2 linclust failed" >> pipeline.log; exit 1; }

# Search with HHblits and HMMER to identify ILP-like sequences
hhblits -i preprocess/ref_clust_rep_seq.fasta -d hhblits_db/ref_db -o preprocess/ref_hhblits.out -n 3 -cpu "$max_cpus" -v 1 || { echo "$(date) - HHblits failed" >> pipeline.log; exit 1; }
hmmsearch -o preprocess/ref_hmm.out --cpu "$max_cpus" input/ilp.hmm preprocess/ref_clust_rep_seq.fasta || { echo "$(date) - HMMER failed" >> pipeline.log; exit 1; }
python combine_hits.py preprocess/ref_hhblits.out preprocess/ref_hmm.out preprocess/ref_candidates_temp.fasta || { echo "$(date) - combine_hits.py failed" >> pipeline.log; exit 1; }

# Deduplicate at 95% identity with CD-HIT
cd-hit -i preprocess/ref_candidates_temp.fasta -o preprocess/ref_candidates.fasta -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date) - CD-HIT failed" >> pipeline.log; exit 1; }

# Build HMM profile and search for additional similarity features
hhmake -i preprocess/ref_candidates.fasta -o preprocess/ref_candidates.hhm -v 0 || { echo "$(date) - hhmake failed" >> pipeline.log; exit 1; }
hhsearch -i preprocess/ref_candidates.hhm -d input/ilp_db.hhm -o preprocess/ref_hhsearch.out -cpu "$max_cpus" -v 1 || { echo "$(date) - HHsearch failed" >> pipeline.log; exit 1; }

# BLAST search against reference ILPs for sequence identity
makeblastdb -in input/ref_ILPs.fasta -dbtype prot -out input/ref_ILPs_db || { echo "$(date) - makeblastdb failed" >> pipeline.log; exit 1; }
blastp -query preprocess/ref_candidates.fasta -db input/ref_ILPs_db -out preprocess/ref_blast.out -outfmt 6 -num_threads "$max_cpus" -evalue 1e-5 || { echo "$(date) - BLAST failed" >> pipeline.log; exit 1; }

# Run InterProScan for domain annotation to assist in cleavage site prediction
interproscan.sh -i preprocess/ref_candidates.fasta -dp -f tsv -iprlookup -goterms -pa -cpu "$max_cpus" > preprocess/ref_interpro.tsv || { echo "$(date) - InterProScan failed" >> pipeline.log; exit 1; }

# Process sequences into prepro, pro, and mature forms using SignalP and InterPro
for seq in $(seqkit seq -n preprocess/ref_candidates.fasta); do
    python preprocess_ilp.py preprocess/ref_candidates.fasta "$seq" prepro "preprocess/prepro_${seq}.fasta" || { echo "$(date) - preprocess_ilp.py failed (prepro)" >> pipeline.log; exit 1; }
    python preprocess_ilp.py preprocess/ref_candidates.fasta "$seq" pro "preprocess/pro_${seq}.fasta" || { echo "$(date) - preprocess_ilp.py failed (pro)" >> pipeline.log; exit 1; }
    python preprocess_ilp.py preprocess/ref_candidates.fasta "$seq" mature "preprocess/mature_${seq}.fasta" || { echo "$(date) - preprocess_ilp.py failed (mature)" >> pipeline.log; exit 1; }
done

# Extract features and generate labels for ML training
python extract_training_features.py preprocess/ref_candidates.fasta preprocess/ref_hhblits.out preprocess/ref_hmm.out preprocess/ref_hhsearch.out preprocess/ref_blast.out preprocess/ref_interpro.tsv preprocess/ > preprocess/ref_features.csv || { echo "$(date) - extract_training_features.py failed" >> pipeline.log; exit 1; }
python generate_labels.py input/ref_ILPs.fasta preprocess/ref_features.csv preprocess/ref_labels.csv || { echo "$(date) - generate_labels.py failed" >> pipeline.log; exit 1; }

# Clean up temporary files with pigz
pigz -f preprocess/ref_clust* preprocess/ref_candidates_temp.fasta
echo "$(date) - 00b_prepare_training.sh completed" >> pipeline.log
touch preprocess/.done
