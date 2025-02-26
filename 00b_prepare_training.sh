#!/bin/bash
# Purpose: Step 00b - Prepare training data for ML by processing reference sequences
# Inputs: input/ref_ILPs.fasta (annotated ILP/non-ILP sequences), input/ilp.hmm, input/ilp_db.hhm
# Outputs: preprocess/ref_features.csv, preprocess/ref_labels.csv, preprocess/*_${type}.fasta
# Config: config.yaml (max_cpus)
# Log: pipeline.log
# Notes: Uses DeepNeuropePred via preprocess_ilp.py for sequence processing
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

max_cpus=$(yq e '.max_cpus' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 00b_prepare_training.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f preprocess/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 00b_prepare_training.sh (already done)" >> pipeline.log
    exit 0
fi

# Check required input files
for file in input/ref_ILPs.fasta input/ilp.hmm input/ilp_db.hhm; do
    if [ ! -f "$file" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: $file not found" >> pipeline.log
        exit 1
    fi
done

mkdir -p preprocess hhblits_db

if [ ! -f hhblits_db/ref_db ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Building HHblits database from UniClust" >> pipeline.log
    hhblits -i input/ref_ILPs.fasta -d uniclav -o /dev/null -cpu "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHblits database build failed" >> pipeline.log; exit 1; }
    mv uniclav* hhblits_db/
    ln -s hhblits_db/uniclav hhblits_db/ref_db
fi

mmseqs linclust input/ref_ILPs.fasta preprocess/ref_clust preprocess/tmp --min-seq-id 0.95 --threads "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 linclust failed" >> pipeline.log; exit 1; }
hhblits -i preprocess/ref_clust_rep_seq.fasta -d hhblits_db/ref_db -o preprocess/ref_hhblits.out -n 3 -cpu "$max_cpus" -v 1 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHblits failed" >> pipeline.log; exit 1; }
hmmsearch -o preprocess/ref_hmm.out --cpu "$max_cpus" input/ilp.hmm preprocess/ref_clust_rep_seq.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HMMER failed" >> pipeline.log; exit 1; }
python combine_hits.py preprocess/ref_hhblits.out preprocess/ref_hmm.out preprocess/ref_candidates_temp.fasta || { echo "$(date '+%Y-%m-%d %H:%M:%S') - combine_hits.py failed" >> pipeline.log; exit 1; }
cd-hit -i preprocess/ref_candidates_temp.fasta -o preprocess/ref_candidates.fasta -c 0.95 -n 5 -T "$max_cpus" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - CD-HIT failed" >> pipeline.log; exit 1; }
hhmake -i preprocess/ref_candidates.fasta -o preprocess/ref_candidates.hhm -v 0 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - hhmake failed" >> pipeline.log; exit 1; }
hhsearch -i preprocess/ref_candidates.hhm -d input/ilp_db.hhm -o preprocess/ref_hhsearch.out -cpu "$max_cpus" -v 1 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - HHsearch failed" >> pipeline.log; exit 1; }
makeblastdb -in input/ref_ILPs.fasta -dbtype prot -out input/ref_ILPs_db || { echo "$(date '+%Y-%m-%d %H:%M:%S') - makeblastdb failed" >> pipeline.log; exit 1; }
blastp -query preprocess/ref_candidates.fasta -db input/ref_ILPs_db -out preprocess/ref_blast.out -outfmt 6 -num_threads "$max_cpus" -evalue 1e-5 || { echo "$(date '+%Y-%m-%d %H:%M:%S') - BLAST failed" >> pipeline.log; exit 1; }
interpro_path=$(yq e '.interpro_path' config.yaml)
"$interpro_path" -i preprocess/ref_candidates.fasta -dp -f tsv -iprlookup -goterms -pa -cpu "$max_cpus" > preprocess/ref_interpro.tsv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - InterProScan failed" >> pipeline.log; exit 1; }

# Process sequences into prepro, pro, mature forms using DeepNeuropePred
seq_count=$(seqkit seq -n preprocess/ref_candidates.fasta | wc -l)
for type in prepro pro mature; do
    for seq in $(seqkit seq -n preprocess/ref_candidates.fasta); do
        python preprocess_ilp.py preprocess/ref_candidates.fasta "$seq" "$type" "preprocess/${type}_${seq}.fasta" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - preprocess_ilp.py failed ($type) for $seq" >> pipeline.log; exit 1; }
    done
    processed_count=$(ls preprocess/${type}_*.fasta | wc -l)
    if [ "$processed_count" -ne "$seq_count" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Incomplete processing for $type (expected $seq_count, got $processed_count)" >> pipeline.log
        exit 1
    fi
done

# Extract sequence-based features
python extract_training_features.py preprocess/ref_candidates.fasta preprocess/ref_hhblits.out preprocess/ref_hmm.out preprocess/ref_hhsearch.out preprocess/ref_blast.out preprocess/ref_interpro.tsv preprocess/ > preprocess/ref_features.csv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - extract_training_features.py failed" >> pipeline.log; exit 1; }
python generate_labels.py input/ref_ILPs.fasta preprocess/ref_features.csv preprocess/ref_labels.csv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_labels.py failed" >> pipeline.log; exit 1; }

pigz -f preprocess/ref_clust* preprocess/ref_candidates_temp.fasta

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 00b_prepare_training.sh completed in ${runtime}s" >> pipeline.log
touch preprocess/.done
