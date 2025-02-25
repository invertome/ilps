# ILP Identification and Phylogenetic Analysis Pipeline

This pipeline identifies insulin-like peptides (ILPs) from transcriptomic data, performs machine learning-based annotation, and conducts phylogenetic analysis across prepropeptide, propeptide, and mature peptide forms. It is optimized for efficiency and scalability, leveraging parallel processing and memory-efficient techniques.

## Prerequisites
- **Tools**: `curl`, `pigz`, `seqkit`, `TransDecoder`, `cd-hit`, `mmseqs2`, `hhblits`, `hmmsearch`, `hhmake`, `hmmbuild`, `hhsearch`, `blastp`, `interproscan.sh`, `signalp6`, `colabfold_batch`, `mafft`, `trimal`, `FastTree`, `iqtree`, `foldtree`, `R` (with `ape` package), `ete3`, `autophy`, `meme`, `ame`, `fimo`, `taxonkit`, `parallel`, `mamba`, `snakemake`, `pymol`, `yq`
- **Python Libraries**: `BioPython`, `pandas`, `scikit-learn`, `xgboost`, `shap`, `matplotlib`, `seaborn`, `logomaker`, `psutil`, `pyyaml`
- **Hardware**: Multi-core CPU recommended; GPU optional for ColabFold

### Installation via Conda
```bash
conda create -n ilp_pipeline python=3.9
conda activate ilp_pipeline
conda install -c bioconda curl pigz seqkit transdecoder cd-hit mmseqs2 hhsuite hmmer blast interproscan mafft trimal fasttree iqtree meme taxonkit parallel mamba snakemake yq
conda install -c conda-forge r-base r-ape ete3 autophy psutil pymol-open-source
conda install -c anaconda pandas scikit-learn xgboost matplotlib seaborn
pip install biopython shap logomaker pyyaml
# ColabFold requires separate installation (see https://github.com/sokrypton/ColabFold)
# FoldTree requires custom installation (see https://github.com/DessimozLab/fold_tree)
```

## Directory Structure
- `input/`: Input FASTA files (e.g., `9606_T1.fasta`)
- `preprocess/`: Preprocessed sequences and reference data
- `candidates/`: Candidate ILP sequences and metadata
- `analysis/`: Phylogenetic trees, structural models, and intermediate files
- `output/`: Final tables and plots (`prepro/`, `pro/`, `mature/` subdirectories)

## Pipeline Steps

### 00a_fetch_references.sh
- **Purpose**: Fetches annotated ILP and non-ILP sequences from UniProt for training.
- **Process**: Queries UniProt for ILPs (e.g., insulin, relaxin) across Metazoa, model organisms, and non-model organisms; fetches an equal number of non-ILPs; annotates with InterProScan and `annotate_references.py`.
- **Output**: `input/ref_ILPs.fasta` with `[ILP]` or `[non-ILP]` tags.
- **Details**: Uses `curl` for API queries, `pigz` for compression, and balances dataset for ML training.

### 00b_prepare_training.sh
- **Purpose**: Prepares reference data for machine learning.
- **Process**: Clusters sequences with `linclust`, searches with `HHblits` and `HMMER`, deduplicates with `CD-HIT`, builds HMM profiles, performs BLAST, and annotates domains with InterProScan. Processes sequences into prepro, pro, and mature forms with `preprocess_ilp.py`. Extracts features (`extract_training_features.py`) and generates labels (`generate_labels.py`).
- **Output**: `preprocess/ref_features.csv`, `preprocess/ref_labels.csv`, `preprocess/*_{type}.fasta`.
- **Details**: Optimizes with multi-threading (`-T $max_cpus`) and chunked processing in Python scripts.

### 01_preprocess.sh
- **Purpose**: Preprocesses input transcriptomes.
- **Process**: Translates nucleotide sequences with `TransDecoder` if needed, filters by length based on reference ILPs (`calc_ref_lengths.py`), deduplicates with `CD-HIT`.
- **Output**: `preprocess/*_preprocessed.fasta`.
- **Details**: Multi-threaded with `seqkit` and `TransDecoder` for efficiency.
  
### 02_run_colabfold.sh
- **Purpose**: Generates structural models for reference and candidate sequences.
- **Process**: Runs ColabFold on prepro, pro, and mature FASTA files from `preprocess/` (training) and `analysis/pdbs/` (candidates), skipping existing PDBs with checkpointing.
- **Inputs**: `preprocess/*_{type}.fasta` (from `00b`), `analysis/pdbs/*_{type}.fasta` (from `05`).
- **Output**: `preprocess/*.pdb`, `analysis/pdbs/*.pdb`.
- **Details**: Uses `parallel` with incremental checks, GPU support optional, validates inputs; runs after training preparation to support feature extraction and again for candidate structures.

### 03_identify_candidates.sh
- **Purpose**: Identifies ILP candidates from preprocessed transcriptomes.
- **Process**: Clusters with `linclust`, searches with `HHblits` and `HMMER`, builds HMM profiles, performs batched BLAST, and annotates domains with InterProScan in parallel.
- **Output**: `candidates/*_candidates.fasta`, metadata files (`*_hhblits.out`, etc.).
- **Details**: Optimizes with `parallel` for InterProScan and batched BLAST for speed.

### 04_annotate_and_novel.sh
- **Purpose**: Annotates candidates with ML predictions and identifies novel ILPs.
- **Process**: Combines candidates into `analysis/all_candidates.fasta`, extracts features with `extract_features.py`, runs ML models (`run_ml.py`) to predict probabilities and novelty.
- **Output**: `analysis/features.csv`, `analysis/predictions.csv`, `analysis/novel_candidates.csv`.
- **Details**: See ML section below for specifics.

### 05_comparative_analysis.sh
- **Purpose**: Performs phylogenetic and structural analysis.
- **Process**: 
  - Filters candidates to ILPs (`filter_ilps.py`).
  - Determines common taxonomy (`determine_common_taxonomy.py`) and filters references (`filter_ref_ilps_by_taxonomy.py`).
  - Processes sequences into prepro, pro, and mature forms (`preprocess_ilp.py`).
  - Aligns with `MAFFT` and trims with `trimal` (cached for efficiency), builds initial trees with `FastTree`, refines with `IQ-TREE`, generates structural trees with `foldtree` in parallel, and creates consensus trees (`consensus_tree.R`).
  - Conducts clade analysis with `ETE` and `Autophy`, followed by motif discovery (`MEME`, `AME`, `FIMO`) and logo generation (`plot_alignment.py`).
- **Output**: `analysis/` (trees, PDBs, clades, plots).
- **Details**: Optimized with `parallel -j 3` to limit concurrent runs, allocating more CPUs per task (e.g., MAFFT, IQ-TREE) for better multi-threading efficiency; caches alignments to speed re-runs.

### 06_generate_outputs.sh
- **Purpose**: Generates final tables and plots.
- **Process**: Runs `generate_tables.py` and `generate_plots.py` for each type (prepro, pro, mature).
- **Output**: `output/*/*.csv` (overview, details, motif enrichment), `output/*/*.png` (counts, heatmap, violin, logos).
- **Details**: Uses taxonkit for taxonomy in plots; chunked processing for memory efficiency.

## Machine Learning Sections (03_annotate_and_novel.sh, run_ml.py)

### Overview
The ML component identifies ILPs and novel candidates using an ensemble of Random Forest (RF) and XGBoost models trained on reference data.

### Feature Extraction (`extract_training_features.py`, `extract_features.py`)
- **Inputs**: Reference (`preprocess/ref_candidates.fasta`) and candidate (`analysis/all_candidates.fasta`) sequences, search outputs (HHblits, HMMER, HHsearch, BLAST), InterPro annotations, and PDBs.
- **Process**: 
  - Extracts sequence similarity scores (HHblits probability, HMMER score, HHsearch probability, BLAST identity).
  - Computes structural similarity (TM-scores) against dynamic and standard references (1TRZ, 6RLX, Bombyxin-II) for prepro, pro, and mature forms using `tmalign`.
  - Adds InterPro domains as categorical features.
  - Uses chunked processing (`chunksize` based on available memory) to handle large datasets efficiently.
- **Output**: Feature matrices (`preprocess/ref_features.csv`, `analysis/features.csv`) with numeric and categorical columns.

### Label Generation (`generate_labels.py`)
- **Inputs**: `input/ref_ILPs.fasta` (annotated references), `preprocess/ref_features.csv`.
- **Process**: Assigns binary labels (1 for ILP, 0 for non-ILP) based on `[ILP]` tags, aligns with feature IDs.
- **Output**: `preprocess/ref_labels.csv`.

### Model Training and Prediction (`run_ml.py`)
- **Inputs**: `analysis/features.csv`, `preprocess/ref_features.csv`, `preprocess/ref_labels.csv`, max CPUs.
- **Process**:
  1. **Data Loading**: Loads features and labels with dynamic chunk sizing based on available memory (`psutil`).
  2. **Initial Training**: Trains a full RF model on all features using `GridSearchCV` for hyperparameter tuning (n_estimators, max_depth, min_samples_split), leveraging all CPUs (`n_jobs=$max_cpus`).
  3. **Feature Selection**: Computes SHAP values to identify the top 10 most important features (e.g., TM-scores, HHblits probability), reducing dimensionality for efficiency and generalization.
  4. **Optimized Training**: Re-trains RF and XGBoost on selected features with tuned hyperparameters (XGBoost: n_estimators, max_depth, learning_rate).
  5. **Prediction**: Combines RF and XGBoost probabilities (average) to predict ILP likelihood; flags novel ILPs (probability > 0.7, HHsearch < 70, BLAST < 30).
  6. **Visualization**: Generates a SHAP summary plot from the full RF model for feature importance insight.
- **Output**: `analysis/predictions.csv` (probabilities), `analysis/novel_candidates.csv` (novelty flags), `output/shap_summary.png`, saved models (`rf_model.joblib`, `xgb_model.joblib`).
- **Details**: Balances accuracy and efficiency with feature selection; uses multi-threading for training; saves models for reuse.

## Usage
1. Place input FASTA files in `input/` (e.g., `9606_T1.fasta`).
2. Run scripts sequentially:
   ```bash
   bash 00a_fetch_references.sh && bash 00b_prepare_training.sh && bash 00c_run_colabfold.sh && bash 01_preprocess.sh && bash 02_identify_candidates.sh && bash 03_annotate_and_novel.sh && bash 04_comparative_analysis.sh && bash 05_generate_outputs.sh
