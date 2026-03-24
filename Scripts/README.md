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

### `shotgun_quick.sh` **[RECOMMENDED FOR METHANE GENE ANALYSIS]**
**Fast methane gene detection pipeline for shotgun metagenomic data**

**Purpose:**
- Rapid detection and quantification of methane metabolism genes in environmental samples
- Identifies methanotrophs (methane oxidizers via pmoA, mmoX) and methanogens (methane producers via mcrA)
- Provides taxonomic context for detected genes
- Generates abundance metrics (raw counts and RPKM-normalized)

**Detailed Workflow:**

```
Raw FASTQ → fastp QC → Host Removal → Kraken2 Taxonomy → DIAMOND Gene Search → Quantification
```

**Step 1: Quality Control with fastp**
- Removes adapters and low-quality bases
- Performs deduplication to reduce PCR bias
- Filters reads shorter than 50bp
- Generates HTML quality reports

**Step 2: Host DNA Removal with Bowtie2**
- Maps reads to plant/animal host genome
- Removes contaminating host DNA
- Retains unmapped microbial reads
- Reduces false-positive gene matches

**Step 3: Taxonomic Profiling with Kraken2**
- Assigns taxonomy using k-mer matching
- Identifies microbial community composition
- Provides context for gene detections
- Generates species-level profiles

**Step 4: Functional Gene Search with DIAMOND**
- Translated BLASTX search against methane gene database
- Identifies: pmoA, pmoB, pmoC (particulate MMO)
           mmoX, mmoY, mmoZ (soluble MMO)
           mcrA, mcrB, mcrG (methyl-coenzyme M reductase)
- Reports best hits with identity >40%, e-value <1e-5
- Extracts gene sequences from reads

**Step 5: Quantification and Normalization**
- Counts reads per gene
- Calculates RPKM (Reads Per Kilobase per Million mapped reads)
- Normalizes for gene length and sequencing depth
- Generates summary tables

**Expected Runtime:**
- Single sample: 50-100 minutes (depends on read count)
- 6 samples in parallel: 60-120 minutes total

**Requirements:**
```bash
# Create environment (one-time setup)
conda create -n quick_search -c bioconda -c conda-forge \
    python=3.10 kraken2 fastp diamond bowtie2

# Activate before running
conda activate quick_search
```

**Required Databases:**
- **Kraken2 Database**
  - Standard: Basic taxonomy (bacteria, archaea, viruses)
  - PlusPF: Includes protozoa and fungi
  - Location: `Data/reference_dbs/kraken2_db/`
  
- **DIAMOND Database**
  - Custom methane gene database
  - Location: `Data/reference_dbs/DIAMOND/methane_master_db.dmnd`

**Input Requirements:**
- Paired-end Illumina reads (FASTQ format, gzipped)
- Naming convention: `{SAMPLE_ID}_R1*.fastq.gz` and `{SAMPLE_ID}_R2*.fastq.gz`
- Location: `Data/raw_data/shotgun/`
- Typical size: 500MB - 2GB per sample (compressed)

**Usage Examples:**

```bash
# Activate environment first
conda activate quick_search

# Single sample - process sample 53394
bash Scripts/01_pipelines/shotgun_quick.sh 53394

# Multiple samples - process in sequence
bash Scripts/01_pipelines/shotgun_quick.sh 53394 53395 53396

# Auto-detect all samples in raw_data/shotgun/
bash Scripts/01_pipelines/shotgun_quick.sh --auto

# Custom thread allocation (default: 4)
bash Scripts/01_pipelines/shotgun_quick.sh --threads 16 53394

# Process with host removal enabled
bash Scripts/01_pipelines/shotgun_quick.sh --threads 8 53394

# Skip specific steps (for troubleshooting)
# Edit script to comment out unwanted steps
```

**Output Files:**

1. **Quality Control Outputs** (`Data/processed_data/fastp_cleaned/`)
   - `{SAMPLE}_R1_clean.fastq.gz` - Quality-filtered forward reads
   - `{SAMPLE}_R2_clean.fastq.gz` - Quality-filtered reverse reads
   - `{SAMPLE}_fastp.html` - Interactive quality report
   - `{SAMPLE}_fastp.json` - Machine-readable QC metrics

2. **Host Removal Outputs** (`Data/processed_data/host_removed/`)
   - `{SAMPLE}_R1_host_removed.fastq.gz` - Host-depleted forward reads
   - `{SAMPLE}_R2_host_removed.fastq.gz` - Host-depleted reverse reads

3. **Taxonomic Profiles** (`Data/processed_data/kraken2_output/`)
   - `{SAMPLE}_kraken2.txt` - Full classification results
   - `{SAMPLE}_kraken2_report.txt` - Taxonomy summary with abundances

4. **Functional Gene Results** (`Data/functional_analysis/methane_genes/`)
   - `{SAMPLE}_combined_methane_hits.txt` - All DIAMOND hits (raw)
   - `{SAMPLE}_filtered_methane_hits.txt` - Filtered hits (identity >40%)
   - `{SAMPLE}_gene_counts.txt` - Read counts per gene
   - `{SAMPLE}_gene_rpkm.txt` - RPKM-normalized abundances
   - `{SAMPLE}_methane_summary.txt` - Summary statistics

5. **Log Files** (`Logs/`)
   - `{SAMPLE}_pipeline.log` - Complete pipeline log
   - `{SAMPLE}_fastp.log` - fastp detailed log
   - `{SAMPLE}_diamond.log` - DIAMOND search log

**Performance Tips:**
- Use `--threads` to match your CPU count
- For 6+ samples, use HPC batch submission (see 02_hpc/)
- Keep compressed FASTQ files to save disk space (~70% reduction)
- Monitor memory usage: ~8-16GB per Kraken2 process

**Troubleshooting:**
- **"Database not found"**: Check paths in script, ensure databases are downloaded
- **"Out of memory"**: Reduce parallel samples or request more RAM on HPC
- **"No hits found"**: Check input file naming, verify read quality
- **"Slow performance"**: Use more threads, check if database is on fast storage

---

### `16S.sh` **[FOR AMPLICON SEQUENCING]**
**QIIME2 pipeline for 16S rRNA amplicon data with checkpoint system**

**Purpose:**
- Process demultiplexed paired-end 16S amplicon sequencing data
- Generate high-quality ASV (Amplicon Sequence Variant) or OTU tables
- Perform taxonomic classification of bacterial and archaeal communities
- Support resumable analysis with checkpoint system

**Detailed Workflow:**

```
Raw FASTQ → Import → QC Visualization → Trim Primers → Denoise/Cluster → Taxonomy → Export
```

**Step 1: Data Import (~2-5 min)**
- Imports paired-end FASTQ files into QIIME2 format (.qza)
- Requires manifest file listing sample IDs and file paths
- Validates file integrity and format
- Creates: `Data/processed_data/demux-paired-end.qza`

**Step 2: Quality Visualization (~1-2 min)**
- Generates interactive quality plots
- Shows per-base quality scores for forward and reverse reads
- Helps determine trimming/truncation parameters
- Creates: `demux-paired-end.qzv` (view at https://view.qiime2.org)

**Step 3: Primer Trimming (~5-10 min, optional)**
- Removes primer sequences from reads
- Common primers: 515F/806R (V4), 341F/805R (V3-V4)
- Uses cutadapt for accurate trimming
- Oriented read trimming (handles mixed orientations)

**Step 4: Denoising with DADA2 (~30-60 min)**
- Removes sequencing errors and chimeras
- Generates exact amplicon sequence variants (ASVs)
- Merges paired-end reads
- Parameters:
  - `--p-trunc-len-f`: Truncate forward reads (e.g., 240)
  - `--p-trunc-len-r`: Truncate reverse reads (e.g., 200)
  - `--p-trim-left-f/r`: Remove primers (if not done in Step 3)
- Creates: `feature-table.qza`, `rep-seqs.qza`, `denoising-stats.qza`

**Alternative: OTU Clustering Mode (~20-40 min)**
- Traditional OTU clustering at 97% identity
- Uses vsearch for clustering
- Faster but less precise than DADA2
- Better for cross-study comparisons

**Step 5: Taxonomic Classification (~10-30 min)**
- Assigns taxonomy to ASVs/OTUs
- Uses pre-trained classifier (Silva, Greengenes, or custom)
- Provides confidence scores for assignments
- Creates: `taxonomy.qza`

**Step 6: Export Results (~2-5 min)**
- Exports QIIME2 artifacts to standard formats
- Feature table → TSV (for R, Python analysis)
- Taxonomy → TSV
- Representative sequences → FASTA

**Checkpoint System:**
- Automatically saves progress after each step
- Resume from last successful step if interrupted
- Skip completed steps on re-run
- Checkpoints stored in: `Logs/checkpoints/`

**Requirements:**
```bash
# Create QIIME2 environment (one-time, ~15 minutes)
conda create -n qiime2 -c qiime2 qiime2

# Activate before running
conda activate qiime2
```

**Required Files:**

1. **Manifest File** (`Data/raw_data/16S/manifest.tsv`)
```tsv
sample-id	forward-absolute-path	reverse-absolute-path
Sample1	$PWD/Data/raw_data/16S/Sample1_R1.fastq.gz	$PWD/Data/raw_data/16S/Sample1_R2.fastq.gz
Sample2	$PWD/Data/raw_data/16S/Sample2_R1.fastq.gz	$PWD/Data/raw_data/16S/Sample2_R2.fastq.gz
```

2. **Taxonomic Classifier** (download pre-trained or train custom)
```bash
# Download Silva 138 classifier (V4 region, 515F/806R)
wget https://data.qiime2.org/2024.10/common/silva-138-99-515-806-nb-classifier.qza \
  -O Data/reference_dbs/silva-138-99-515-806-nb-classifier.qza
```

**Usage Examples:**

```bash
# Activate QIIME2 environment
conda activate qiime2

# Interactive mode (prompts at each step)
bash Scripts/01_pipelines/16S.sh

# Non-interactive/batch mode for HPC
bash Scripts/01_pipelines/16S.sh --non-interactive --mode denoise

# Use OTU clustering instead of DADA2
bash Scripts/01_pipelines/16S.sh --mode cluster

# Resume from checkpoint after interruption
bash Scripts/01_pipelines/16S.sh
# (automatically detects completed steps)
```

**Expected Runtime:**
- Small dataset (10 samples, 50K reads each): 30-60 minutes
- Medium dataset (50 samples, 100K reads each): 1-3 hours
- Large dataset (200 samples, 200K reads each): 4-8 hours

**Output Files:**

1. **QIIME2 Artifacts** (`Data/processed_data/`)
   - `demux-paired-end.qza` - Imported sequences
   - `feature-table.qza` - ASV/OTU abundance table
   - `rep-seqs.qza` - Representative sequences
   - `taxonomy.qza` - Taxonomic assignments
   - `*.qzv` - Visualization files (drag to view.qiime2.org)

2. **Exported Tables** (`Results/denoise_mode/` or `Results/cluster_mode/`)
   - `exported-table/feature-table.tsv` - Read counts per ASV per sample
   - `exported-taxonomy/taxonomy.tsv` - Taxonomy per ASV
   - `exported-seqs/dna-sequences.fasta` - Representative sequences

3. **Statistics** 
   - `denoising-stats.tsv` - Reads retained at each filtering step
   - Typical retention: 70-90% of input reads

**Performance Tips:**
- **Memory**: 4-8GB per sample during DADA2
- **Storage**: ~1-2GB per 100 samples (artifacts)
- **CPUs**: DADA2 uses 1 thread per sample (batch parallelization)
- Use `--non-interactive` mode for HPC job submission

**Troubleshooting:**
- **"Low read retention after DADA2"**: 
  - Check quality plots, adjust truncation lengths
  - Ensure primers are properly trimmed
  - Verify reads overlap (sum of trunc lengths < amplicon length)
  
- **"Taxonomic classifier not found"**:
  - Download appropriate classifier for your primer set
  - Place in `Data/reference_dbs/`
  
- **"Out of memory"**:
  - Process smaller batches
  - Request more RAM on HPC (8GB minimum, 16GB recommended)
  
- **"Checkpoint error"**:
  - Clear checkpoints: `rm -rf Logs/checkpoints/*`
  - Re-run from beginning

**Key Parameters to Adjust:**

For **DADA2 mode**:
```bash
# In the script, modify these based on quality plots:
--p-trunc-len-f 240    # Truncate forward reads (keep high-quality region)
--p-trunc-len-r 200    # Truncate reverse reads (usually lower quality)
--p-trim-left-f 0      # Remove primers if not already trimmed
--p-trim-left-r 0      # 
--p-max-ee-f 2         # Maximum expected errors (forward)
--p-max-ee-r 2         # Maximum expected errors (reverse)
```

For **OTU clustering mode**:
```bash
--p-perc-identity 0.97  # 97% identity threshold (standard)
--p-threads 4           # Parallel threads for clustering
```

---

### `shotgun_MAGs.sh` **[FOR METAGENOME-ASSEMBLED GENOMES]**
**QIIME2-moshpit shotgun metagenomic analysis pipeline for MAG reconstruction**

**Purpose:**
- Comprehensive analysis of shotgun metagenomic data
- Reconstruct Metagenome-Assembled Genomes (MAGs) from community DNA
- Perform assembly, binning, and functional annotation
- Identify novel microbial genomes from environmental samples

**Detailed Workflow:**

```
Raw FASTQ → QC → Assembly → Gene Prediction → Binning → MAG Refinement → Annotation
```

**Step 1: Data Import (~5-10 min)**
- Imports paired-end shotgun reads into QIIME2
- Requires manifest file with sample paths
- Validates file format and integrity

**Step 2: Quality Control with FastQC + fastp (~20-40 min)**
- **FastQC** generates detailed quality reports (before/after)
- **fastp** performs:
  - Adapter trimming (auto-detection)
  - Quality filtering (Q20 threshold)
  - Length filtering (>50bp)
  - Deduplication
- **MultiQC** aggregates all QC reports into single dashboard
- Typical retention: 80-95% of reads

**Step 3: Metagenomic Assembly (~2-12 hours)**

**Option A: MEGAHIT** (default)
- Memory-efficient de Bruijn graph assembler
- Suitable for complex metagenomes
- Uses iterative k-mer strategy (21, 29, 39, 59, 79, 99, 119, 141)
- Memory: ~5-20GB per sample
- Speed: ~2-4 hours for 10M read pairs

**Option B: SPAdes (metaSPAdes)**
- More accurate, better for low-abundance organisms
- Requires more memory (~40-100GB)
- Slower than MEGAHIT (~4-8 hours)
- Better contiguity (longer contigs)

**Assembly Modes:**
- **per-sample**: Assemble each sample independently
  - Preserves strain-level diversity
  - More MAGs but potentially fragmented
  
- **coassembly**: Combine all samples before assembly
  - Better for low-coverage organisms
  - Produces consensus genomes
  - Loses strain-level resolution

**Step 4: Gene Prediction with Prodigal (~30-60 min)**
- Identifies protein-coding genes in contigs
- Uses metagenomics mode (`-p meta`)
- Outputs: gene sequences (nucleotide and amino acid)

**Step 5: Taxonomic Classification (~1-2 hours)**
- Classifies contigs using CAT/BAT or Kaiju
- Assigns taxonomy to assembled sequences
- Provides confidence scores

**Step 6: Functional Annotation (~2-4 hours)**
- **DIAMOND** against KEGG/COG/Pfam databases
- Identifies metabolic pathways
- Annotates: 
  - Carbohydrate metabolism
  - Nitrogen/sulfur cycling
  - Energy metabolism
  - Antibiotic resistance genes

**Step 7: Binning and MAG Reconstruction (~1-3 hours)**
- **MetaBAT2** or **CONCOCT** for automated binning
- Groups contigs into genome bins based on:
  - Tetranucleotide frequency
  - Coverage profiles across samples
  - Co-abundance patterns
- Generates draft MAG genomes

**Step 8: MAG Quality Assessment (~30 min)**
- **CheckM** evaluates completeness and contamination
- High-quality MAG: >90% complete, <5% contamination
- Medium-quality: >50% complete, <10% contamination
- Reports: rRNA genes, tRNA counts, strain heterogeneity

**Requirements:**
```bash
# Create specialized environment (one-time, ~30 minutes)
conda create -n qiime2-moshpit -c qiime2 qiime2-moshpit megahit metabat2 checkm-genome

# Activate before running
conda activate qiime2-moshpit
```

**Database Requirements:**
- **KEGG**: Functional annotation (~10GB)
- **CheckM**: Quality assessment (~1.5GB)
- **CAT/BAT**: Taxonomic classification (~25GB)

**Configuration Options:**
Edit script to set:
```bash
ASSEMBLER="megahit"              # or "spades"
ASSEMBLY_MODE="per-sample"       # or "coassembly"
MIN_CONTIG_LENGTH=1000          # Filter short contigs
BINNING_MIN_COVERAGE=2          # Minimum contig coverage
```

**Usage:**
```bash
# Activate environment
conda activate qiime2-moshpit

# Run full pipeline
bash Scripts/01_pipelines/shotgun_MAGs.sh

# Monitor progress (run in another terminal)
tail -f Logs/shotgun_pipeline.log
```

**Expected Runtime:**
- Single sample (10M reads): 6-12 hours
- Five samples (coassembly): 12-24 hours
- Large dataset (100M+ reads): 24-48 hours

**Expected Resource Usage:**
- **MEGAHIT**: 8-20GB RAM, 4-8 CPUs
- **SPAdes**: 40-100GB RAM, 16-32 CPUs
- **Disk space**: 50-200GB per sample (includes all intermediates)

**Output Files:**

1. **Assembly** (`Results/assembly/`)
   - `final_contigs.fasta` - Assembled contigs (>1000bp)
   - `assembly_stats.txt` - N50, L50, total length
   - Typical metrics: N50 5-50 Kb, longest contig 100-500 Kb

2. **Gene Predictions** (`Results/gene_prediction/`)
   - `predicted_genes.fna` - Gene sequences (nucleotide)
   - `predicted_genes.faa` - Protein sequences (amino acid)
   - `gene_annotations.gff` - Gene coordinates

3. **Functional Annotations** (`Results/functional_annotation/`)
   - `kegg_annotations.txt` - KEGG pathway assignments
   - `cog_annotations.txt` - COG functional categories
   - `pathway_abundance.txt` - Metabolic pathway completeness

4. **MAG Bins** (`Results/binning/`)
   - `bin.001.fa` to `bin.XXX.fa` - Individual MAG genomes
   - `bins_summary.txt` - Completeness/contamination per bin
   - Typical output: 5-50 MAGs per sample

5. **Quality Reports** (`Results/qc/`)
   - `multiqc_report.html` - Comprehensive QC dashboard
   - `checkm_results.txt` - MAG quality metrics

**Interpretation Guide:**

**Assembly Quality Metrics:**
- Good: N50 > 10 Kb, max contig > 100 Kb
- Poor: N50 < 2 Kb, fragmented assembly
- **Fix**: Increase read depth, try coassembly

**MAG Quality Standards (MIMAG):**
- **High-quality**: >90% complete, <5% contamination, 23S/16S/5S rRNA, >18 tRNAs
- **Medium-quality**: ≥50% complete, <10% contamination
- **Low-quality**: <50% complete, <10% contamination

**Performance Tips:**
- Start with MEGAHIT for initial exploration
- Use SPAdes for high-value samples requiring best quality
- Coassembly works best with 3-10 related samples
- Keep intermediate files for troubleshooting (consumes disk space)
- Use HPC for datasets >20M read pairs

**Troubleshooting:**
- **"Assembly failed / Out of memory"**:
  - Use MEGAHIT instead of SPAdes
  - Reduce k-mer sizes
  - Subsample reads for testing
  
- **"Few/no bins generated"**:
  - Check assembly quality (N50 too low)
  - Ensure multiple samples for coverage profiles
  - Adjust binning parameters (coverage thresholds)
  
- **"Low MAG completeness"**:
  - Common for complex communities
  - Try coassembly with more samples
  - Increase sequencing depth
  
- **"High contamination in MAGs"**:
  - Strain heterogeneity (multiple strains)
  - Run MAG refinement (manual curation)
  - Use DAS Tool for bin aggregation

---

## 02_downstream/ - Downstream Analysis

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
bash Scripts/01_pipelines/shotgun_quick.sh 53394 53395 53396

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
bash Scripts/01_pipelines/16S.sh --non-interactive

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
│   │   ├── 16S/                 # 16S amplicon FASTQ files
│   │   └── shotgun/              # Shotgun metagenomic FASTQ files
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
