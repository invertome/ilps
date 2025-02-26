# ILP Identification and Phylogenetic Analysis Pipeline

This pipeline identifies insulin-like peptides (ILPs) from transcriptomic data, performs machine learning-based annotation, and conducts phylogenetic analysis across prepropeptide, propeptide, and mature peptide forms. It is optimized for efficiency and scalability, leveraging parallel processing and memory-efficient techniques.

## Prerequisites
- **Tools**: `curl`, `pigz`, `seqkit`, `TransDecoder`, `cd-hit`, `mmseqs2`, `hhblits`, `hmmsearch`, `hhmake`, `hmmbuild`, `hhsearch`, `blastp`, `interproscan.sh`, `signalp6`, `colabfold_batch`, `mafft`, `trimal`, `FastTree`, `iqtree`, `foldtree`, `R` (with `ape` package), `ete3`, `autophy`, `meme`, `ame`, `fimo`, `taxonkit`, `parallel`, `mamba`, `snakemake`, `pymol`, `yq`
- **Python Libraries**: `BioPython`, `pandas`, `scikit-learn`, `xgboost`, `shap`, `matplotlib`, `seaborn`, `logomaker`, `psutil`, `pyyaml`
- **Hardware**: Multi-core CPU recommended; GPU optional for ColabFold

### Installation via Conda

```bash
# Create and activate the base environment with Python 3.11 for SignalP 6
conda create -n ilp_pipeline python=3.11
conda activate ilp_pipeline

# Add necessary channels
conda config --append channels bioconda
conda config --append channels conda-forge

# Install core bioinformatics tools
conda install -c bioconda curl pigz seqkit transdecoder cd-hit mmseqs2 hhsuite blast interproscan mafft trimal fasttree iqtree meme taxonkit parallel mamba snakemake yq

# Install additional dependencies, including SignalP 6 requirements
conda install -c conda-forge biopython r-base r-ape ete3 psutil pymol-open-source pillow numpy matplotlib tqdm pytorch=1.13
conda install -c anaconda pandas scikit-learn xgboost seaborn

# Install Python libraries via pip
pip install shap logomaker pyyaml
pip install 'autophy @ git+https://github.com/aortizsax/autophy@main'

# Install SignalP 6 manually (requires license from DTU Health Tech)
# Download signalp-6.0h.fast.tar.gz from DTU Health Tech, then:
# tar xfv signalp-6.0h.fast.tar.gz
# cd signalp-6-package
# pip install .
# ColabFold requires separate installation (see https://github.com/sokrypton/ColabFold)
# FoldTree is integrated via Snakemake in the pipeline (see https://github.com/DessimozLab/fold_tree)
```

## Directory Structure
- `input/`: Input FASTA files (e.g., `9606_T1.fasta`)
- `preprocess/`: Preprocessed sequences and reference data
- `candidates/`: Candidate ILP sequences and metadata
- `analysis/`: Phylogenetic trees, structural models, and intermediate files
- `output/`: Final tables and plots (`prepro/`, `pro/`, `mature/` subdirectories)

## Pipeline Steps

### 00a_fetch_references.sh

- **Purpose**: Fetches annotated ILP and non-ILP sequences from UniProt and generates HMM profiles for training and candidate identification.
- **Process**: Queries UniProt for ILPs (e.g., insulin, relaxin) across Metazoa and curated references (e.g., P01308), balances with non-ILPs, annotates with InterProScan and `annotate_references.py`, aligns ILP sequences with MAFFT, builds `ilp.hmm` with `hmmbuild`, and `ilp_db.hhm` with `hhmake`.
- **Output**: `input/ref_ILPs.fasta` (annotated), `input/ilp.hmm`, `input/ilp_db.hhm`.
- **Details**: Uses `curl` for API calls, `pigz` for compression, multi-threaded alignment; includes dependency checks.

### 00b_prepare_training.sh

- **Purpose**: Prepares reference data for machine learning with sequence-based features.
- **Process**: Clusters sequences with `linclust`, searches with `HHblits` and `HMMER`, deduplicates with `CD-HIT`, performs BLAST, and annotates domains with InterProScan. Processes sequences into prepro, pro, and mature forms with `preprocess_ilp.py`. Extracts sequence features (`extract_training_features.py`) and generates labels (`generate_labels.py`).
- **Inputs**: `input/ref_ILPs.fasta`, `input/ilp.hmm`, `input/ilp_db.hhm`.
- **Output**: `preprocess/ref_features.csv`, `preprocess/ref_labels.csv`, `preprocess/*_{type}.fasta`.
- **Details**: Multi-threaded (`-T $max_cpus`), chunked processing, validates inputs and sequence completion; structural features deferred to candidate processing.

### 01_preprocess.sh

- **Purpose**: Preprocesses input transcriptomes.
- **Process**: Translates nucleotide sequences with `TransDecoder` if needed, filters by length based on reference ILPs (`calc_ref_lengths.py`), deduplicates with `CD-HIT`.
- **Inputs**: `input/[0-9]*_*.fasta`, `input/ref_ILPs.fasta`.
- **Output**: `preprocess/*_preprocessed.fasta`.
- **Details**: Multi-threaded with `seqkit` and `TransDecoder`, validates input files.

### 02_identify_candidates.sh

- **Purpose**: Identifies ILP candidates from preprocessed transcriptomes.
- **Process**: Clusters with `linclust`, searches with `HHblits` and `HMMER`, builds HMM profiles, performs batched BLAST, and annotates domains with InterProScan in parallel.
- **Inputs**: `preprocess/[0-9]*_preprocessed.fasta`, `input/ref_ILPs.fasta`, `input/ilp.hmm`, `input/ilp_db.hhm`.
- **Output**: `candidates/*_candidates.fasta`, metadata files (`*_hhblits.out`, etc.).
- **Details**: Optimizes with `parallel`, caches BLAST/InterProScan, validates inputs.

### 03_annotate_and_novel.sh

- **Purpose**: Performs initial ML annotation of candidates based on sequence features.
- **Process**: Combines candidates into `analysis/all_candidates.fasta`, extracts sequence-based features with `extract_features.py`, runs ML models (`run_ml.py`) for initial probabilities and novelty prediction.
- **Inputs**: `candidates/[0-9]*_candidates.fasta`, `preprocess/ref_features.csv`, `preprocess/ref_labels.csv`.
- **Output**: `analysis/all_candidates.fasta`, `analysis/predictions.csv`, `analysis/novel_candidates.csv`.
- **Details**: Initial pass without structural features; see ML section for specifics.

### 04_run_colabfold.sh

- **Purpose**: Generates structural models for identified ILP candidates only.
- **Process**: Runs ColabFold on prepro, pro, and mature forms of candidate sequences from `candidates/*_candidates.fasta`, processed via `preprocess_ilp.py`, skipping existing PDBs with checkpointing.
- **Inputs**: `candidates/*_candidates.fasta`.
- **Output**: `analysis/pdbs/*.pdb`.
- **Details**: Uses `parallel` with incremental checks, GPU support optional, validates inputs; predicts structures only for candidates to optimize computational efficiency.

### 05_comparative_analysis.sh

- **Purpose**: Performs phylogenetic and structural analysis with comprehensive consensus trees.
- **Process**: 
  - Filters candidates to ILPs (`filter_ilps.py`), determines taxonomy (`determine_common_taxonomy.py`), filters references (`filter_ref_ilps_by_taxonomy.py`).
  - Aligns with `MAFFT`, trims with `trimal`, builds sequence trees with `FastTree` and `IQ-TREE`, generates structural trees with `foldtree` (Foldtree, LDDT, TM metrics) using candidate PDBs.
  - Creates consensus trees (`consensus_tree_with_support.R`) for all combinations:
    - Per type: Sequence + Foldtree, Sequence + LDDT, Sequence + TM, Foldtree + LDDT, Foldtree + TM, LDDT + TM, Foldtree + LDDT + TM, Sequence + Foldtree + LDDT + TM.
    - Across types: Sequence and Foldtree combinations (e.g., prepro + pro).
  - Conducts clade analysis with `ETE` and `Autophy`, motif discovery (`MEME`, `AME`, `FIMO`), and logo generation (`plot_alignment.py`).
- **Inputs**: `analysis/all_candidates.fasta`, `analysis/predictions.csv`, `analysis/novel_candidates.csv`, `input/[0-9]*_*.fasta`, `input/ref_ILPs.fasta`, `analysis/pdbs/*.pdb`.
- **Output**: `analysis/` (trees like `*_consensus_seq_foldtree.tre`, PDBs, clades, plots).
- **Details**: Uses `parallel -j 3`, caches alignments, integrates Foldtree via Snakemake, validates `mamba`/`snakemake`.

### 06_generate_outputs.sh

- **Purpose**: Generates final outputs with structural features for manuscript preparation.
- **Process**: 
  - Re-runs `extract_features.py` with structural features from candidate PDBs, generates tables (`generate_tables.py`) and plots (`generate_plots.py`) for each type.
  - Produces annotated FASTA files and metadata TSV (`generate_output_fasta_and_metadata.py`).
  - Creates 3D structure figures (`generate_structure_figures.py`) for candidates.
- **Inputs**: `analysis/all_candidates.fasta`, `analysis/predictions.csv`, `analysis/novel_candidates.csv`, `candidates/*_blast.out`, `candidates/*_interpro.tsv`, `clades_ete_{type}/`, `clades_autophy_{type}/`, `analysis/ilp_candidates.fasta`, `analysis/ref_ILPs_filtered.fasta`, `input/[0-9]*_*.fasta`, `analysis/pdbs/`.
- **Output**: `output/*/*.csv` (overview, details, motif enrichment), `output/*/*.png` (counts, heatmap, violin, logos), `output/*/ilps.fasta` (annotated ILPs), `output/comparative_metadata.tsv` (sequence metadata), `output/figures/*.png` (3D structures).
- **Details**: Uses `taxonkit` for taxonomy, validates PyMOL, chunked processing; includes final feature extraction with structural data.

## Machine Learning Sections (03_annotate_and_novel.sh, 06_generate_outputs.sh, run_ml.py)

### Overview

Identifies ILPs and novel candidates using an ensemble of Random Forest (RF) and XGBoost models trained on reference data, with an initial sequence-based pass and a final pass incorporating structural features.

### Feature Extraction (`extract_training_features.py`, `extract_features.py`)

- **Initial Pass (03_annotate_and_novel.sh)**:
  - **Inputs**: Reference (`preprocess/ref_candidates.fasta`) and candidate (`analysis/all_candidates.fasta`) sequences, search outputs (HHblits, HMMER, HHsearch, BLAST), InterPro annotations.
  - **Process**: 
    - Extracts sequence similarity scores (HHblits probability, HMMER score, HHsearch probability, BLAST identity).
    - Adds physicochemical properties (hydrophobicity, charge) and InterPro domains.
    - Uses chunked processing based on available memory (`psutil`).
  - **Output**: `preprocess/ref_features.csv`, `analysis/features_initial.csv`.
- **Final Pass (06_generate_outputs.sh)**:
  - **Inputs**: `analysis/all_candidates.fasta`, `analysis/pdbs/*.pdb`.
  - **Process**: 
    - Repeats sequence feature extraction as above.
    - Computes structural similarity (TM-scores, pLDDT) against dynamic and standard references (1TRZ, 6RLX, Bombyxin-II) for prepro, pro, and mature forms using `tmalign` on candidate PDBs.
  - **Output**: `analysis/features_final.csv`.

### Label Generation (`generate_labels.py`)

- **Inputs**: `input/ref_ILPs.fasta`, `preprocess/ref_features.csv`.
- **Process**: Assigns binary labels (1 for ILP, 0 for non-ILP) based on `[ILP]` tags.
- **Output**: `preprocess/ref_labels.csv`.

### Model Training and Prediction (`run_ml.py`)

- **Inputs**: `analysis/features_initial.csv`, `preprocess/ref_features.csv`, `preprocess/ref_labels.csv`, max CPUs.
- **Process**:
  1. Loads data with dynamic chunk sizing (`psutil`).
  2. Trains RF with `GridSearchCV` for hyperparameter tuning, computes SHAP values for feature selection (top 10 features).
  3. Re-trains RF and XGBoost with 5-fold cross-validation, reporting AUC, precision, recall (`output/ml_metrics.txt`).
  4. Predicts ILP probabilities (RF + XGBoost average), flags novel ILPs (probability > 0.7, HHsearch < 70, BLAST < 30).
  5. Generates SHAP summary plot.
- **Output**: `analysis/predictions.csv`, `analysis/novel_candidates.csv`, `output/shap_summary.png`, `output/rf_model.joblib`, `output/xgb_model.joblib`, `output/ml_metrics.txt`.
- **Details**: Multi-threaded training, reusable models, robust validation; initial pass in `03` uses sequence features only, with structural refinement possible post-`04` if re-run.


## Parameter Customization
- **Thresholds**: Key thresholds such as `ilp_prob_threshold`, `hhsearch_novel_threshold`, and `blast_novel_threshold` can be adjusted in `config.yaml` to fine-tune ILP identification and novelty detection.
- **Tools and Paths**: Ensure all tool paths (e.g., `interpro_path`, `colabfold_path`) are correctly set in `config.yaml`.

## Error Handling
- The pipeline includes checks for missing tools, input files, and failed commands. Logs are recorded in `pipeline.log` for troubleshooting.
- For external tool failures (e.g., `Foldtree` in `05_comparative_analysis.sh`), check the corresponding log files in `analysis/` or `fold_tree/` directories.

  
## Usage Example
1. Place your transcriptome FASTA files in `input/` with TaxID prefixes (e.g., `9606_T1.fasta`).
2. Update `config.yaml` with appropriate tool paths and parameters.
3. Run the pipeline sequentially:
   ```bash
   bash 00a_fetch_references.sh && bash 00b_prepare_training.sh && bash 01_preprocess.sh && bash 02_identify_candidates.sh && bash 03_annotate_and_novel.sh && bash 04_run_colabfold.sh && bash 05_comparative_analysis.sh && bash 06_generate_outputs.sh


   
