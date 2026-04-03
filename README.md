# Shotgun Metagenomic Analysis Pipeline

Fast methane metabolism gene detection pipeline for shotgun metagenomic data.

**Focus**: Detection and quantification of methane metabolism genes (pmoA, mmoX, mcrA) in environmental samples.

---

## Quick Start

```bash
# 1. Create conda environment
conda env create -f quick_search.yml
conda activate quick_search

# 2. Run pipeline (auto-detect all samples)
bash Scripts/shotgun_quick.sh --auto --threads 16

# 3. Or submit to HPC
sbatch Scripts/run_shotgun.slurm
---

## Pipeline Overview

**Workflow:**
```
Raw FASTQ → fastp QC → Host Removal → Kraken2 Taxonomy → DIAMOND Gene Search → Quantification
**Steps:**
1. **fastp QC** - Quality control, adapter trimming, deduplication
2. **Bowtie2 Host Removal** - Remove plant/animal DNA contamination (optional)
3. **Kraken2 Taxonomy** - Taxonomic classification of microbial community
4. **DIAMOND** - Functional gene search against methane gene database
5. **Quantification** - Count genes and calculate RPKM-normalized abundance

**Expected Runtime:**
- Single sample (~10M reads): 15-20 minutes
- 6 samples (sequential): 90-120 minutes

---

## Installation

### 1. Create Conda Environment

```bash
# Create environment from YAML file
conda env create -f quick_search.yml

# Activate environment
conda activate quick_search
```

**Installed packages:**
- Python 3.10
- Quality Control: fastp, bbmap
- Alignment: bowtie2
- Taxonomy: kraken2
- Functional Analysis: diamond, hmmer, orfm, singlem, smafa, pplacer

### 2. Prepare Databases

**Kraken2 Database** (required)
```bash
# Option 1: Download pre-built standard database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
mkdir -p Data/reference_dbs
tar -xzf k2_standard_*.tar.gz -C Data/reference_dbs/

# Option 2: Build custom database
kraken2-build --standard --threads 16 --db Data/reference_dbs
**DIAMOND Database** (required)
- Place methane gene database at: `Data/reference_dbs/DIAMOND/methane_master_db.dmnd`

**Bowtie2 Host Index** (optional)
- Place host genome index at: `Data/reference_dbs/host_genomes/plant_host.*.bt2`
- If not present, pipeline will skip host removal step

---

## Usage

### Command Line

```bash
# Activate environment first
conda activate quick_search

# Auto-detect all samples in Data/raw_data/
bash Scripts/shotgun_quick.sh --auto --threads 16

# Process specific samples
bash Scripts/shotgun_quick.sh 53394 53395 53396 --threads 16

# Read sample IDs from file
bash Scripts/shotgun_quick.sh --sample-list samples.txt --threads 16

# Force re-run all steps (skip checkpoints)
bash Scripts/shotgun_quick.sh --auto --force --threads 16

# View help
bash Scripts/shotgun_quick.sh --help
### HPC Batch Job (SLURM)

```bash
# Submit job
sbatch Scripts/run_shotgun.slurm

# Monitor job
squeue -u $USER

# Check logs
tail -f Logs/slurm_<JOB_ID>.out
cat Logs/slurm_<JOB_ID>.err

# Check job status
sacct -j <JOB_ID> --format=JobID,State,ExitCode,Elapsed,MaxRSS
**Default HPC settings** (in run_shotgun.slurm):
- Partition: large_336
- Nodes: 1
- CPUs: 16
- Memory: 128GB
- Time limit: 48 hours

---

## Input Data

### File Requirements

- **Format**: Paired-end Illumina reads (FASTQ, gzipped)
- **Location**: `Data/raw_data/`
- **Naming**: `{SAMPLE_ID}_R1*.fastq.gz` and `{SAMPLE_ID}_R2*.fastq.gz`
  - Examples: 
    - `53394_R1_A15.fastq.gz`, `53394_R2_A15.fastq.gz`
    - `53395_R1_A16.fastq.gz`, `53395_R2_A16.fastq.gz`
- **Size**: Typically 500MB - 2GB per sample (compressed)

### Sample Detection

The `--auto` flag will automatically find all samples matching the pattern `*_R1*.fastq.gz` in `Data/raw_data/`.

---

## Output Files

All outputs organized by sample ID:

### 1. Quality Control (`Data/processed_data/fastp_cleaned/`)
- `{SAMPLE}_R1_cleaned.fastq.gz` - Quality-filtered forward reads
- `{SAMPLE}_R2_cleaned.fastq.gz` - Quality-filtered reverse reads
- `{SAMPLE}_fastp.html` - Interactive QC report
- `{SAMPLE}_fastp.json` - Machine-readable QC metrics

### 2. Host Removal (`Data/processed_data/host_removed/`)
- `{SAMPLE}_R1_host_removed.fastq.gz` - Host-depleted forward reads
- `{SAMPLE}_R2_host_removed.fastq.gz` - Host-depleted reverse reads
- `{SAMPLE}_host_removal_stats.txt` - Bowtie2 alignment statistics

### 3. Taxonomic Classification (`Results/kraken2_output/`)
- `kraken2_output_{SAMPLE}/kraken2_output.txt` - Per-read classifications
- `kraken2_output_{SAMPLE}/kraken2_report.txt` - Taxonomic summary with abundances

### 4. Functional Gene Analysis (`Results/functional_analysis/`)
- `{SAMPLE}_combined_methane_hits.txt` - All DIAMOND hits (raw)
- `{SAMPLE}_filtered_methane_hits.txt` - Filtered methane genes only
- `{SAMPLE}_gene_counts.txt` - Read counts per gene
- `{SAMPLE}_gene_rpkm.txt` - RPKM-normalized gene abundances
- `{SAMPLE}_methane_summary.txt` - Summary statistics

### 5. Logs (`Logs/`)
- `{SAMPLE}_quick_search.log` - Complete pipeline log per sample
- `{SAMPLE}_fastp.log` - fastp detailed output
- `slurm_{JOB_ID}.out` - SLURM job standard output
- `slurm_{JOB_ID}.err` - SLURM job error output

---

## Results Interpretation

### Quality Control Metrics

**Good quality indicators:**
- Q20 bases: >95%
- Q30 bases: >85%
- Reads passing filter: >90%
- Duplication rate: <5%

### Taxonomic Classification

**Typical results:**
- Classification rate: 20-40% (environmental samples)
- Unclassified: 60-80% (normal for complex metagenomes)

**Methanotrophic bacteria to look for:**
- **Type I** (Gammaproteobacteria): Methylomonas, Methylococcus, Methylomicrobium
- **Type II** (Alphaproteobacteria): Methylocystis, Methylosinus
- **Type X**: Methylocaldum

### Functional Gene Hits

**DIAMOND output interpretation:**
- **Total hits**: All alignments to methane gene database
- **Filtered hits**: After removing ribosomal/non-target genes
- **MMO**: Methane monooxygenase (soluble form)
- **pMMO**: Particulate methane monooxygenase
- **MCR**: Methyl-coenzyme M reductase (methanogens)

**RPKM values:**
- Normalizes for gene length and sequencing depth
- Comparable across samples
- Higher RPKM = higher gene abundance

---

## Directory Structure

```
shotgun/
├── quick_search.yml              # Conda environment specification
├── README.md                     # This file
│
├── Scripts/
│   ├── shotgun_quick.sh          # Main pipeline script
│   └── run_shotgun.slurm         # SLURM job submission script
│
├── Data/
│   ├── raw_data/                 # Input FASTQ files
│   │   ├── 53394_R1_A15.fastq.gz
│   │   ├── 53394_R2_A15.fastq.gz
│   │   └── ...
│   │
│   ├── processed_data/
│   │   ├── fastp_cleaned/        # Quality-controlled reads
│   │   └── host_removed/         # Host-depleted reads
│   │
│   ├── reference_dbs/
│   │   ├── hash.k2d              # Kraken2 database files
│   │   ├── taxo.k2d              # (can also be in kraken2_db/)
│   │   ├── opts.k2d
│   │   │
│   │   ├── DIAMOND/
│   │   │   └── methane_master_db.dmnd
│   │   │
│   │   └── host_genomes/         # Optional
│   │       └── plant_host.*.bt2
│   │
│   ├── temp/                     # Temporary files
│   └── metadata/
│       └── metadata.tsv          # Sample metadata
│
├── Results/                      # Analysis outputs
│   ├── kraken2_output/           # Taxonomic classification results
│   │   └── kraken2_output_{SAMPLE}/
│   └── functional_analysis/      # DIAMOND gene search results
│       ├── {SAMPLE}_combined_methane_hits.txt
│       ├── {SAMPLE}_filtered_methane_hits.txt
│       ├── {SAMPLE}_gene_counts.txt
│       ├── {SAMPLE}_gene_rpkm.txt
│       └── {SAMPLE}_methane_summary.txt
│
└── Logs/                         # All log files
    ├── slurm_*.out
    ├── slurm_*.err
    └── *_quick_search.log
```
---
## Troubleshooting

### Common Issues

**"No samples detected with --auto"**
- Check file naming: must match `*_R1*.fastq.gz` pattern
- Verify files are in `Data/raw_data/` directory
- Try: `ls Data/raw_data/*_R1*.fastq.gz`

**"Database not found"**
```bash
# Check Kraken2 database
ls Data/reference_dbs/*.k2d
# or
ls Data/reference_dbs/kraken2_db/*.k2d

# Check DIAMOND database
ls Data/reference_dbs/DIAMOND/methane_master_db.dmnd
**"Tool not found" or "Command not found"**
```bash
# Ensure environment is activated
conda activate quick_search

# Verify tools are installed
conda list | grep -E "fastp|diamond|kraken2|bowtie2"

# Recreate environment if tools are missing
conda env remove -n quick_search
conda env create -f quick_search.yml
**"Out of memory" on HPC**
- Increase memory in `run_shotgun.slurm`: `#SBATCH --mem=256G`
- Reduce parallel samples: `MAX_PARALLEL_SAMPLES=1`
- Check Kraken2 memory usage (typically needs 100GB+ for standard DB)

**"Job terminated without completing all samples"**
```bash
# Check job status
sacct -j <JOB_ID> --format=JobID,State,ExitCode,Elapsed,MaxRSS

# Common causes:
# - Time limit exceeded: increase --time
# - Memory exceeded: increase --mem
# - Node failure: resubmit job
**"Host removal skipped"**
- This is normal if Bowtie2 index not found
- Pipeline will continue without host removal
- To enable: create index at `Data/reference_dbs/host_genomes/plant_host`

**"No methane genes detected"**
- Verify DIAMOND database is correct
- Check sample quality (fastp HTML report)
- View DIAMOND output: `cat Results/functional_analysis/{SAMPLE}_combined_methane_hits.txt`
- Sample may genuinely lack methanotrophs/methanogens

---

## Performance Optimization

### Threading
```bash
# Adjust based on available CPUs
bash Scripts/shotgun_quick.sh --auto --threads 32
### Parallel Processing
Edit `Scripts/shotgun_quick.sh` or set environment variable:
```bash
export MAX_PARALLEL_SAMPLES=4  # Process 4 samples simultaneously
```
Requires GNU parallel: `conda install -c conda-forge parallel`

### Disk Space Management
```bash
# Remove intermediate files after successful run
rm -rf Data/temp/*
rm -rf Data/processed_data/fastp_cleaned/  # Keep only final results
### Checkpoints
Pipeline automatically skips completed steps. To force re-run:
```bash
bash Scripts/shotgun_quick.sh --auto --force
---

## Database Information

### Kraken2 Database Options

**Standard Database** (~45GB download, ~100GB RAM needed)
- Bacteria, archaea, viruses, and human
- Good for general microbial communities

**PlusPF Database** (~75GB download, ~150GB RAM needed)
- Adds protozoa and fungi
- Better for complex environmental samples

**MiniKraken** (~8GB download, ~8GB RAM needed)
- Smaller, faster, less accurate
- Good for testing or resource-limited systems

### DIAMOND Database

Custom methane metabolism gene database should contain:
- **Methanotroph genes**: pmoA, pmoB, pmoC, mmoX, mmoY, mmoZ
- **Methanogen genes**: mcrA, mcrB, mcrG
- Curated from NCBI, UniProt, or custom sources

---

## Citations

If you use this pipeline, please cite:

- **fastp**: Chen et al. (2018) Bioinformatics
- **Kraken2**: Wood et al. (2019) Genome Biology  
- **DIAMOND**: Buchfink et al. (2021) Nature Methods
- **Bowtie2**: Langmead & Salzberg (2012) Nature Methods

---

## Contact

**Author**: Jiayi Chen  
**Email**: jc224@imperial.ac.uk  
**Last Updated**: April 2026  
**Version**: 1.0

For issues or questions, please contact the author.

---

## Changelog

### Version 1.0 (April 2026)
- Initial release
- Supports auto-detection of samples
- Integrated QC, taxonomy, and functional analysis
- SLURM batch job support
- Checkpoint system for resumable analysis
