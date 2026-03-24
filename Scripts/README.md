# Metagenomic Analysis Pipeline Scripts

Scripts for shotgun metagenomic analysis with focus on methane metabolism genes (pmoA, mmoX, mcrA) and microbial community profiling.

---

## Directory Structure

```
Scripts/
├── 01_pipelines/         # Main analysis pipelines
├── 02_hpc/              # HPC job submission scripts
└── 03_downstream/       # Downstream data analysis and visualization
```

---

## 01_pipelines/ - Main Analysis Pipelines

### `quick_search.sh` **[RECOMMENDED]**
**Fast methane gene detection pipeline for shotgun metagenomic data**

**Purpose:**
- Rapid detection and quantification of methane metabolism genes
- Identifies methanotrophs (methane oxidizers) and methanogens (methane producers)

**Workflow:**
1. **fastp** - Quality control and adapter trimming
2. **Bowtie2** - Host DNA removal (optional)
3. **Kraken2** - Taxonomic profiling
4. **DIAMOND** - Functional gene search (pmoA, mmoX, mcrA, etc.)
5. **Quantification** - Gene abundance and RPKM normalization

**Requirements:**
```bash
conda create -n quick_search -c bioconda -c conda-forge \
    python=3.10 kraken2 fastp diamond bowtie2
```

**Databases:**
- Kraken2: `Data/reference_dbs/kraken2_db/` (standard or PlusPF)
- DIAMOND: `Data/reference_dbs/DIAMOND/methane_master_db.dmnd`

**Usage:**
```bash
# Activate environment
conda activate quick_search

# Single sample
bash Scripts/01_pipelines/quick_search.sh 53394

# Multiple samples
bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396

# Auto-detect all samples
bash Scripts/01_pipelines/quick_search.sh --auto

# Custom threads
bash Scripts/01_pipelines/quick_search.sh --threads 16 53394
```

**Outputs:**
- `Data/processed_data/fastp_cleaned/` - Quality-controlled reads
- `Data/processed_data/host_removed/` - Host-removed reads
- `Data/processed_data/kraken2_output/` - Taxonomic profiles
- `Data/functional_analysis/methane_genes/` - Gene search results

---

### `main.sh`
**QIIME2 pipeline for amplicon data with checkpoint system**

**Purpose:**
- Process demultiplexed paired-end FASTQ files
- Generate OTU/ASV tables and representative sequences

**Features:**
- Checkpoint system for resumable analysis
- Supports DADA2 denoising or OTU clustering
- Interactive and non-interactive (batch) modes

**Modes:**
- `denoise` - DADA2 denoising (default)
- `cluster` - OTU clustering

**Requirements:**
```bash
conda create -n qiime2 -c qiime2 qiime2
```

**Usage:**
```bash
# Interactive mode
bash Scripts/01_pipelines/main.sh

# Non-interactive (HPC batch mode)
bash Scripts/01_pipelines/main.sh --non-interactive --mode denoise
```

---

### `shotgun.sh`
**QIIME2 shotgun metagenomic analysis pipeline**

**Purpose:**
- Comprehensive shotgun data analysis using QIIME2-moshpit
- Assembly, binning, and functional annotation

**Workflow:**
1. Data import
2. Quality control (FastQC + fastp + MultiQC)
3. Assembly (MEGAHIT or SPAdes)
4. Functional annotation

**Requirements:**
```bash
conda create -n qiime2-moshpit -c qiime2 qiime2-moshpit
```

**Assembly Options:**
- `megahit` - Fast, memory-efficient (default)
- `spades` - More accurate, higher memory

**Assembly Modes:**
- `per-sample` - Individual sample assemblies
- `coassembly` - Combined assembly

**Usage:**
```bash
conda activate qiime2-moshpit
bash Scripts/01_pipelines/shotgun.sh
```

---

## 02_hpc/ - HPC Job Submission Scripts

### PBS Scripts (Imperial HPC - CX3)

#### `run_quick_search.pbs`
Submit quick_search pipeline to PBS scheduler

```bash
qsub Scripts/02_hpc/run_quick_search.pbs
```

**Resources:**
- Nodes: 1
- CPUs: 24
- Memory: 128GB
- Walltime: 24 hours

**Configuration:**
- 5 samples in parallel
- 4 threads per sample (5 × 4 = 20 cores)

---

#### `run_diamond_only.pbs`
Run DIAMOND functional gene search only (skip QC and Kraken2)

```bash
# Default sample
qsub Scripts/02_hpc/run_diamond_only.pbs

# Custom sample
qsub -v SAMPLE_ID=53395 Scripts/02_hpc/run_diamond_only.pbs
```

**Resources:**
- CPUs: 8
- Memory: 32GB
- Walltime: 6 hours

---

#### `run_singlem_only.pbs`
Run SingleM taxonomic profiling only

```bash
qsub -v SAMPLE_ID=53394 Scripts/02_hpc/run_singlem_only.pbs
```

**Resources:**
- CPUs: 8
- Memory: 32GB
- Walltime: 12 hours

---

#### `kraken2_job.sh`
Batch Kraken2 taxonomic profiling for all host-removed samples

```bash
qsub Scripts/02_hpc/kraken2_job.sh
```

**Resources:**
- CPUs: 24
- Memory: 150GB
- Walltime: 24 hours

---

#### `build_plant_index.sh`
Build Bowtie2 index for plant host genome removal

```bash
qsub Scripts/02_hpc/build_plant_index.sh
```

**Input:** `combined_host.fasta`
**Output:** `plant_host.*.bt2` index files

**Resources:**
- CPUs: 8
- Memory: 64GB
- Walltime: 2 hours

---

#### `submit_mothur.sh`
Submit mothur analysis for pmoA gene analysis

```bash
qsub Scripts/02_hpc/submit_mothur.sh
```

**Requirements:**
- Conda environment: `bio_env`
- Mothur batch file: `../01_pipelines/run_pmoA.mth`

---

### SLURM Scripts (General HPC Systems)

#### `submit_quick_search.sbatch`
SLURM submission script for quick_search pipeline

```bash
# Single sample
sbatch submit_quick_search.sbatch 53394

# Multiple samples
sbatch submit_quick_search.sbatch 53394 53395 53396

# Auto-detect
sbatch submit_quick_search.sbatch --auto
```

**Resources:**
- CPUs: 16
- Memory: 64GB
- Time: 24 hours
- Partition: large_336

---

#### `submit_shotgun.sbatch`
SLURM submission for shotgun metagenomic analysis

```bash
sbatch submit_shotgun.sbatch
```

**Resources:**
- CPUs: 16
- Memory: 64GB
- Time: 48 hours

---

#### `submit_qiime.sbatch`
SLURM submission for QIIME2 amplicon analysis

```bash
sbatch submit_qiime.sbatch
```

**Resources:**
- CPUs: 8
- Memory: 64GB
- Time: 24 hours

---

### Utility Scripts

#### `hpc_upload.sh`
Upload pipeline to HPC cluster and auto-submit Slurm job

```bash
bash hpc_upload.sh user@hpc.edu /scratch/user/project
```

**Features:**
- Rsync files to HPC
- Auto-generates Slurm submission script
- Configurable resources

---

## 03_downstream/ - Downstream Analysis

### Python Scripts

#### `analyze_methane_genes.py`
Analyze DIAMOND methane gene search results

**Features:**
- Parse DIAMOND blastx output
- Extract gene and species information
- Classify genes into functional categories (pMMO, sMMO, MCR, etc.)
- Calculate gene abundance and RPKM

**Usage:**
```bash
python analyze_methane_genes.py sample_id
```

**Outputs:**
- Gene-species mapping
- Functional category counts
- Abundance tables

---

#### `visualize_results.py`
Create publication-quality visualizations

**Plots Generated:**
1. Gene category contribution (stacked bar chart)
2. Species contribution (DIAMOND vs SingleM comparison)
3. Comparative analyses

**Features:**
- High-resolution plots (300 DPI)
- Publication-ready figures
- Comparative analysis support

**Usage:**
```bash
python visualize_results.py sample_id output_dir
```

---

### R Scripts

#### `table_combine.R`
Combine QIIME2 feature table with taxonomy annotations

**Input:**
- `Results/denoise_mode/exported-table/feature-table.tsv`
- `Results/denoise_mode/exported-taxonomy/taxonomy.tsv`

**Output:**
- `Results/denoise_mode/taxonomy-abundance-table.tsv`
- `Results/denoise_mode/final-table-with-ranks.tsv`

**Features:**
- Parse QIIME2 taxonomy format
- Merge abundance and taxonomy
- Add taxonomic rank columns

---

#### `taxonomic_analysis.R`
Comprehensive taxonomic analysis and visualization

**Analyses:**
1. Community composition heatmap
2. PCoA ordination
3. PERMANOVA (statistical testing)
4. Alpha and beta diversity

**Features:**
- Filters Cyanobacteria and Chloroplast
- Supports multiple metadata groups
- Publication-quality plots

**Outputs:**
- Heatmaps
- PCoA plots
- Distance matrices
- Statistical test results

---

#### `mob_analysis.R`
Methanotrophic bacteria (MOB) specific analysis

**Target Genera:**
- Type I: Methylobacter, Methylomonas, Methylomicrobium, Methylosoma, Methylosarcina, Methylovulum, Methylogaea, Methylocaldum, Methylococcus
- Type II: Methylocystis, Methylosinus, Methylocapsa, Methylocella
- Type X: Methylocaldum

**Analyses:**
1. MOB identification and abundance
2. Type I vs Type II comparison
3. Taxonomic distribution
4. Statistical testing (Kruskal-Wallis)

**Outputs:**
- MOB abundance tables
- Stacked bar plots
- Box plots
- Statistical summaries

**Datasets Supported:**
- Bacteria only
- Bacteria + Archaea

---

#### `ggpicrust.R`
PICRUSt2 functional pathway analysis

**Purpose:**
- Extract and visualize methane-related pathways
- Identify key enzymes in methane metabolism

**Target Pathways:**
- Methanogenesis: METHANOGENESIS-PWY, PWY-7084, etc.
- Methanotrophy: METHANE-OXIDATION-PWY, PWY-7357, etc.

**Target Enzymes:**
- **Methanogenesis:**
  - MCR (EC 2.8.4.1): Methyl-coenzyme M reductase
  - KO: K00399, K00401, K00577, etc.
- **Methanotrophy:**
  - MMO (EC 1.14.13.25, 1.14.18.3): Methane monooxygenase
  - KO: K10944-K10946, K16157-K16159

**Features:**
- Auto-download enzyme descriptions
- Filter by pathways or enzymes
- Generate abundance tables
- Publication-ready plots

---

## Common Workflows

### Workflow 1: Quick Methane Gene Search
```bash
# 1. Activate environment
conda activate quick_search

# 2. Run pipeline
bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396

# 3. Analyze results
python Scripts/03_downstream/analyze_methane_genes.py 53394

# 4. Visualize
python Scripts/03_downstream/visualize_results.py 53394 Results/plots/
```

### Workflow 2: HPC Batch Processing
```bash
# Submit job
qsub Scripts/02_hpc/run_quick_search.pbs

# Monitor
qstat -u $USER

# Check results
ls Data/functional_analysis/methane_genes/
```

### Workflow 3: QIIME2 Amplicon Analysis
```bash
# 1. Run main pipeline
conda activate qiime2
bash Scripts/01_pipelines/main.sh --non-interactive

# 2. Export tables
qiime tools export ...

# 3. Downstream analysis
Rscript Scripts/03_downstream/table_combine.R
Rscript Scripts/03_downstream/taxonomic_analysis.R
Rscript Scripts/03_downstream/mob_analysis.R
```

---

## Required Conda Environments

### quick_search
```bash
conda create -n quick_search -c bioconda -c conda-forge \
    python=3.10 kraken2 fastp diamond bowtie2
```

### qiime2
```bash
conda create -n qiime2 -c qiime2 qiime2
```

### qiime2-moshpit
```bash
conda create -n qiime2-moshpit -c qiime2 qiime2-moshpit
```

### bio_env (for mothur)
```bash
conda create -n bio_env -c bioconda mothur
```

---

## Expected Directory Structure

```
meta/
├── Data/
│   ├── raw_data/
│   │   └── shotgun/              # Raw FASTQ files
│   ├── processed_data/
│   │   ├── fastp_cleaned/        # QC'd reads
│   │   ├── host_removed/         # Host-removed reads
│   │   └── kraken2_output/       # Taxonomic profiles
│   ├── functional_analysis/
│   │   └── methane_genes/        # DIAMOND results
│   ├── reference_dbs/
│   │   ├── kraken2_db/           # Kraken2 database
│   │   ├── DIAMOND/              # DIAMOND databases
│   │   └── host_genomes/         # Host reference genomes
│   └── metadata/
│       └── metadata.tsv          # Sample metadata
├── Results/                       # Analysis outputs
├── Logs/                         # Log files
└── Scripts/                      # This directory
```

---

## Database Setup

### Kraken2 Database
```bash
# Option 1: Build standard database (28GB, slow)
kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db

# Option 2: Download pre-built (faster)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
mkdir -p Data/reference_dbs/kraken2_db
tar -xzf k2_standard_*.tar.gz -C Data/reference_dbs/kraken2_db/
```

### DIAMOND Database
Place custom methane gene database at:
`Data/reference_dbs/DIAMOND/methane_master_db.dmnd`

---

## Author & Version

**Author:** Jiayi Chen (jc224@imperial.ac.uk)  
**Last Updated:** March 2026  
**Pipeline Version:** 1.0  

---

## References

- **QIIME2:** Bolyen et al. (2019) Nature Biotechnology
- **Kraken2:** Wood et al. (2019) Genome Biology
- **DIAMOND:** Buchfink et al. (2015) Nature Methods
- **fastp:** Chen et al. (2018) Bioinformatics

---

## Notes

1. Always activate the appropriate conda environment before running scripts
2. Ensure databases are downloaded and properly configured
3. Check resource requirements before submitting HPC jobs
4. Monitor job logs in `Logs/` directory for errors
5. Adjust thread counts based on available resources

---

## Troubleshooting

**Issue:** Conda environment not found
```bash
# Verify environment exists
conda env list

# Recreate if needed
conda create -n quick_search ...
```

**Issue:** Database not found
```bash
# Check database paths
ls -lh Data/reference_dbs/kraken2_db/
ls -lh Data/reference_dbs/DIAMOND/

# Verify database integrity
kraken2 --db Data/reference_dbs/kraken2_db/ --inspect
```

**Issue:** Out of memory on HPC
- Increase memory allocation in PBS/SLURM scripts
- Reduce number of parallel samples
- Use swap space if available

---

For questions or issues, contact: jc224@imperial.ac.uk
