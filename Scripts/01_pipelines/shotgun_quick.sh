#!/bin/bash

################################################################################
# Quick Metagenomic Search Pipeline - Methane Gene Detection
################################################################################
# 
# PURPOSE:
#   Fast detection and quantification of methane metabolism genes in shotgun
#   metagenomic data. Identifies methanotrophs (methane oxidizers) and 
#   methanogens (methane producers).
#
# WORKFLOW:
#   1. fastp - Quality control and adapter trimming
#   2. Host Removal - Remove host contamination (optional, if reference available)
#   3. Kraken2 - Taxonomic profiling using 16S marker genes
#   4. DIAMOND - Functional gene search (pmoA, mmoX, mcrA, etc.)
#   5. Quantification - Gene abundance and RPKM normalization
#
# REQUIREMENTS:
#   • Conda environment: quick_search
#     - fastp (fast quality control)
#     - kraken2 (taxonomic classification)
#     - DIAMOND (functional gene search)
#     - bowtie2 (host removal)
#   
#   • Databases (auto-downloaded or manual):
#     - Kraken2 database: Data/reference_dbs/kraken2_db/ (Standard or PlusPF)
#     - DIAMOND database: Data/reference_dbs/DIAMOND/methane_master_db.dmnd
#
#   • Input data:
#     - Paired-end shotgun reads: Data/raw_data/shotgun/<SAMPLE>_R1*.fastq.gz
#                                                        <SAMPLE>_R2*.fastq.gz
#
# SETUP:
#   1. Create conda environment:
#      conda create -n quick_search -c bioconda -c conda-forge \
#          python=3.10 kraken2 fastp diamond bowtie2
#   
#   2. Download Kraken2 database (8-28GB depending on version):
#      conda activate quick_search
#      kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db
#      # Or download pre-built: wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_*.tar.gz
#   
#   3. Copy DIAMOND database to Data/reference_dbs/DIAMOND/
#
# USAGE EXAMPLES:
#   # Activate environment first
#   conda activate quick_search
#
#   # Single sample
#   bash Scripts/01_pipelines/quick_search.sh 53394
#
#   # Multiple samples
#   bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396
#
#   # From sample list file
#   bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt
#
#   # Auto-detect all samples
#   bash Scripts/01_pipelines/quick_search.sh --auto
#
#   # Custom threads (default: 4)
#   bash Scripts/01_pipelines/quick_search.sh --threads 16 53394
#
# HPC USAGE:
#   sbatch Scripts/02_hpc/submit_quick_search.sbatch 53394
#
# OUTPUT:
#   Data/processed_data/fastp_cleaned/         - Quality-controlled reads
#   Data/processed_data/host_removed/          - Host-decontaminated reads (if enabled)
#   Data/processed_data/kraken2_output/        - Taxonomic profiles
#   Data/functional_analysis/methane_genes/    - Gene search results
#   Logs/                                      - Sample-specific logs
#
# AUTHOR: Jiayi
# DATE: 2026-02-01
################################################################################


# ============================================================================
# Environment Check and Configuration
# ============================================================================

# Set temporary directory to a writable location to avoid issues on HPC systems with /tmp
TMPDIR_BASE="${BASE_DIR:-$(pwd)}/Data/temp"
mkdir -p "$TMPDIR_BASE"
export TMPDIR="$TMPDIR_BASE"
export TEMP="$TMPDIR_BASE"
export TMP="$TMPDIR_BASE"
# ======================================================================


# Default configuration
THREADS=${THREADS:-4}              # Number of threads (can be set via environment variable)
MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-1}  # Number of samples to process in parallel
RAW_DATA_SUBDIR="${RAW_DATA_SUBDIR:-shotgun}"  # Subdirectory for raw data (shotgun or empty for root)

# Initialize conda (simplified - HPC scripts should handle this)
# For local use, ensure conda is initialized in your shell
if ! command -v conda &> /dev/null; then
    # Minimal initialization attempt for local use
    [ -f "${HOME}/.bashrc" ] && source "${HOME}/.bashrc" 2>/dev/null
    CONDA_BASE="${CONDA_BASE:-$(conda info --base 2>/dev/null)}"
    [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ] && source "${CONDA_BASE}/etc/profile.d/conda.sh" 2>/dev/null
fi

# Check if conda environment is activated; if running under HPC scheduler, auto-activate
if [[ "$CONDA_DEFAULT_ENV" != "quick_search" ]]; then
    if [[ -n "$SLURM_JOB_ID" ]] || [[ -n "$PBS_JOBID" ]] || [[ -n "$PBS_JOBNAME" ]]; then
        echo "HPC batch job detected. Auto-activating quick_search environment..."
        conda activate quick_search
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to auto-activate quick_search environment"
            exit 1
        fi
    else
        echo "========================================================================"
        echo "WARNING: Conda environment 'quick_search' is not activated"
        echo "========================================================================"
        echo ""
        echo "Please activate the environment first:"
        echo "  conda activate quick_search"
        echo ""
        echo "If you don't have the environment yet, create it:"
        echo "  conda create -n quick_search -c bioconda -c conda-forge \\"
        echo "      python=3.10 fastp diamond kraken2 bowtie2"
        echo ""
        exit 1
    fi
fi

# Auto-detect base directory if not set
if [ -z "$BASE_DIR" ]; then
    # Try to find the project root by looking for Scripts directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
fi

# Parse command line arguments
SAMPLE_LIST=()
SAMPLE_FILE=""
AUTO_DETECT=false

while [[ $# -gt 0 ]]; do
    case $1 in
        step0|--setup)
            echo "Setup Instructions:"
            echo "1. conda create -n quick_search -c bioconda -c conda-forge python=3.10 kraken2 fastp diamond bowtie2"
            echo "2. conda activate quick_search"
            echo "3. Download Kraken2 DB: kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db"
            echo "4. Place DIAMOND database at: Data/reference_dbs/DIAMOND/methane_master_db.dmnd"
            exit 0
            ;;
        --sample-list)
            SAMPLE_FILE="$2"
            shift 2
            ;;
        --auto)
            AUTO_DETECT=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --parallel)
            MAX_PARALLEL_SAMPLES="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --help|-h)
            cat << 'EOF'
========================================================================
  Quick Metagenomic Search - Methane Gene Detection Pipeline
========================================================================

USAGE:
  bash Scripts/01_pipelines/quick_search.sh [OPTIONS] [SAMPLE_IDs...]

OPTIONS:

    step0, --setup        Show recommended conda/package setup (custom env name allowed)
    --sample-list FILE    Read sample IDs from file (one per line)
    --auto                Auto-detect all samples in raw_data directory
    --threads N           Number of threads (default: 4)
    --base-dir PATH       Base directory path (default: auto-detect)
    --help, -h            Show this help message

EXAMPLES:
  # Single sample
  bash Scripts/01_pipelines/quick_search.sh 53394

  # Multiple samples
  bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396

  # From sample list file
  bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt

  # Auto-detect all samples
  bash Scripts/01_pipelines/quick_search.sh --auto

  # Custom thread count
  bash Scripts/01_pipelines/quick_search.sh --threads 16 53394

  # HPC submission
  sbatch Scripts/02_hpc/submit_quick_search.sbatch 53394

REQUIREMENTS:
  1. Conda environment: quick_search (must be activated)
  2. Input data: Data/raw_data/shotgun/<SAMPLE>_R1*.fastq.gz
  3. Databases:
     - Kraken2: Data/reference_dbs/kraken2_db/ (contains hash.k2d, opts.k2d, taxo.k2d)
     - DIAMOND: Data/reference_dbs/DIAMOND/methane_master_db.dmnd

SETUP (first time):
  # Create environment
  conda create -n quick_search -c bioconda -c conda-forge \\
      python=3.10 kraken2 fastp diamond

  # Activate environment
  conda activate quick_search

  # Download Kraken2 database
  kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db
  # Or: wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_*.tar.gz && tar -xzf k2_standard_*.tar.gz -C Data/reference_dbs/kraken2_db/

OUTPUT:
  Data/processed_data/fastp_cleaned/      - Quality-controlled reads
  Data/processed_data/host_removed/       - Host-decontaminated reads (if enabled)
  Data/processed_data/kraken2_output/     - Taxonomic profiles
  Data/functional_analysis/methane_genes/ - Gene detection results
  Logs/                                   - Execution logs

For more information, see: QUICK_SEARCH_USAGE.md
========================================================================
EOF
            exit 0
            ;;
        -*)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
        *)
            SAMPLE_LIST+=("$1")
            shift
            ;;
    esac
done

# Determine samples to process
if [ "$AUTO_DETECT" = true ]; then
    echo "Auto-detecting samples from Data/raw_data/${RAW_DATA_SUBDIR}..."
    RAW_DATA_DIR="${BASE_DIR}/Data/raw_data/${RAW_DATA_SUBDIR}"
    if [ -d "$RAW_DATA_DIR" ]; then
        # Find all unique sample IDs by looking for _R1.fastq.gz files
        while IFS= read -r file; do
            SAMPLE_ID=$(basename "$file" | sed 's/_R1.fastq.gz//')
            SAMPLE_LIST+=("$SAMPLE_ID")
        done < <(find "$RAW_DATA_DIR" -name "*_R1.fastq.gz" -type f)
    fi
elif [ -n "$SAMPLE_FILE" ]; then
    echo "Reading samples from sample list file..."
    if [ ! -f "$SAMPLE_FILE" ]; then
        echo "ERROR: Sample file not found"
        echo "  Specified file: $(basename "$SAMPLE_FILE")"
        echo "  Please check:"
        echo "    1. File path is correct"
        echo "    2. File exists and is readable"
        exit 1
    fi
    while IFS= read -r line || [ -n "$line" ]; do
        # Skip empty lines and comments
        [[ -z "$line" || "$line" =~ ^#.*$ ]] && continue
        SAMPLE_LIST+=("$line")
    done < "$SAMPLE_FILE"
fi

# Check if any samples provided
if [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "========================================================================"
    echo "ERROR: No samples provided"
    echo "======================================================================="
    echo ""
    echo "Usage: bash Scripts/01_pipelines/quick_search.sh [OPTIONS] [SAMPLE_IDs...]"
    echo ""
    echo "Examples:"
    echo "  bash Scripts/01_pipelines/quick_search.sh 53394"
    echo "  bash Scripts/01_pipelines/quick_search.sh --auto"
    echo "  bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt"
    echo ""
    echo "Use --help for more information"
    exit 1
fi

echo "========================================================================"
echo "Quick Metagenomic Search Pipeline"
echo "========================================================================"
echo "Samples to process: ${#SAMPLE_LIST[@]}"
echo "Threads: ${THREADS}"
echo "Base directory: ${BASE_DIR}"
echo "Conda environment: ${CONDA_DEFAULT_ENV:-not detected}"
echo "======================================================================="
echo ""

# ============================================================================
# Prerequisites Check
# ============================================================================

echo "Checking prerequisites for quick_search environment..."
echo ""

# Check required tools in quick_search environment
# Note: kraken2 will be checked later when we switch to kraken_env
MISSING_TOOLS=()

if ! command -v fastp &> /dev/null; then
    MISSING_TOOLS+=("fastp")
fi

if ! command -v diamond &> /dev/null; then
    MISSING_TOOLS+=("diamond")
fi

if ! command -v orfm &> /dev/null; then
    MISSING_TOOLS+=("orfm (orfM)")
fi

if ! command -v hmmsearch &> /dev/null; then
    MISSING_TOOLS+=("hmmsearch (HMMER)")
fi

if ! command -v bc &> /dev/null; then
    MISSING_TOOLS+=("bc")
fi

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo "========================================================================"
    echo "ERROR: Missing required tools in quick_search environment"
    echo "========================================================================"
    echo ""
    echo "The following tools are not found:"
    for tool in "${MISSING_TOOLS[@]}"; do
        echo "  - $tool"
    done
    echo ""
    echo "Please install them in the quick_search environment:"
    echo "  conda activate quick_search"
    echo "  conda install -c bioconda -c conda-forge fastp diamond orfm hmmer"
    echo ""
    echo "Note: kraken2 should be installed separately in kraken_env:"
    echo "  conda create -n kraken_env -c bioconda kraken2"
    echo ""
    exit 1
fi

echo "  ✓ fastp found: $(fastp --version 2>&1 | head -1)"
echo "  ✓ diamond found: $(diamond --version 2>&1 | head -1)"
echo "  ✓ orfm found: $(which orfm)"
echo "  ✓ hmmsearch found: $(which hmmsearch)"
echo "Note: kraken2 will be checked when switching to kraken_env"
echo ""

# Set paths
RAW_DATA_DIR="${BASE_DIR}/Data/raw_data/${RAW_DATA_SUBDIR}"
FASTP_DIR="${BASE_DIR}/Data/processed_data/fastp_cleaned"
HOST_REMOVED_DIR="${BASE_DIR}/Data/processed_data/host_removed"
KRAKEN2_DIR="${BASE_DIR}/Data/processed_data/kraken2_output"
DIAMOND_DIR="${BASE_DIR}/Data/functional_analysis/methane_genes"
LOG_DIR="${BASE_DIR}/Logs"

# Create output directories
mkdir -p ${FASTP_DIR}
mkdir -p ${HOST_REMOVED_DIR}
mkdir -p ${KRAKEN2_DIR}
mkdir -p ${DIAMOND_DIR}
mkdir -p ${LOG_DIR}


# Set database paths - auto-detect or use environment variables
# Note: fastp auto-detects adapters, no need for adapter file

# Auto-detect Kraken2 DB: check kraken2_db/ subdirectory first, fallback to reference_dbs/ itself
if [ -n "$KRAKEN2_DB" ]; then
    :  # user already set it
elif [ -f "${BASE_DIR}/Data/reference_dbs/kraken2_db/hash.k2d" ]; then
    KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs/kraken2_db"
elif [ -f "${BASE_DIR}/Data/reference_dbs/hash.k2d" ]; then
    KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs"
else
    KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs/kraken2_db"  # default for error message
fi
METHANE_DB="${METHANE_DB:-${BASE_DIR}/Data/reference_dbs/DIAMOND/methane_master_db}"

# Host genome reference (if available)
HOST_GENOME_DIR="${BASE_DIR}/Data/reference_dbs/host_genomes"

# Check databases
echo "Checking databases..."
echo ""

DB_MISSING=false

# Check Kraken2 database
if [ -d "$KRAKEN2_DB" ] && [ -f "${KRAKEN2_DB}/hash.k2d" ]; then
    echo "  ✓ Kraken2 database found: ${KRAKEN2_DB}"
else
    echo "  ❌ Kraken2 database NOT found"
    echo "     Expected directory: $KRAKEN2_DB"
    echo "     Must contain: hash.k2d, opts.k2d, taxo.k2d files"
    DB_MISSING=true
fi

# Check DIAMOND database
if [ -f "${METHANE_DB}.dmnd" ]; then
    echo "  ✓ DIAMOND database found: ${METHANE_DB}.dmnd"
elif [ -f "${METHANE_DB}" ]; then
    echo "  ✓ DIAMOND database found: ${METHANE_DB}"
else
    echo "  ❌ DIAMOND database NOT found"
    echo "     Expected: ${METHANE_DB}.dmnd"
    DB_MISSING=true
fi

if [ "$DB_MISSING" = true ]; then
    echo ""
    echo "========================================================================"
    echo "ERROR: Required databases are missing"
    echo "======================================================================="
    echo ""
    echo "To download Kraken2 database (8-28GB):"
    echo "  conda activate quick_search"
    echo "  # Option 1: Build standard database"
    echo "  kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db"
    echo ""
    echo "  # Option 2: Download pre-built database (faster)"
    echo "  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz"
    echo "  mkdir -p Data/reference_dbs/kraken2_db"
    echo "  tar -xzf k2_standard_*.tar.gz -C Data/reference_dbs/kraken2_db/"
    echo ""
    echo "For DIAMOND database:"
    echo "  Copy or create methane_master_db.dmnd in Data/reference_dbs/DIAMOND/"
    echo ""
    exit 1
fi

echo ""
echo "All prerequisites checked - ready to start!"
echo ""

# ============================================================================
# Sample Processing
# ============================================================================

################################################################################
# Function: Process single sample
################################################################################
process_sample() {
    local SAMPLE_ID=$1
    local START_TIME=$(date +%s)
    local LOG_FILE
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_quick_search.log"

    (
    echo "========================================================================"
    echo "Processing Sample: ${SAMPLE_ID}"
    echo "========================================================================"
    echo "Start time: $(date)"
    echo "======================================================================="
    
    # Input files - support both naming patterns
    # Pattern 1: SAMPLE_R1.fastq.gz (standard)
    # Pattern 2: SAMPLE_R1_LABEL.fastq.gz (with sample label)
    if [ -f "${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz" ]; then
        R1_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz"
        R2_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    else
        # Try to find files with pattern SAMPLE_R1_*.fastq.gz
        R1_RAW=$(find "${RAW_DATA_DIR}" -name "${SAMPLE_ID}_R1_*.fastq.gz" -o -name "${SAMPLE_ID}_R1.fastq.gz" | head -1)
        R2_RAW=$(find "${RAW_DATA_DIR}" -name "${SAMPLE_ID}_R2_*.fastq.gz" -o -name "${SAMPLE_ID}_R2.fastq.gz" | head -1)
    fi

    # fastp output files
    R1_CLEAN="${FASTP_DIR}/${SAMPLE_ID}_R1_cleaned.fastq.gz"
    R2_CLEAN="${FASTP_DIR}/${SAMPLE_ID}_R2_cleaned.fastq.gz"
    FASTP_JSON="${FASTP_DIR}/${SAMPLE_ID}_fastp.json"
    FASTP_HTML="${FASTP_DIR}/${SAMPLE_ID}_fastp.html"

    # Host removal output files
    R1_HOST_REMOVED="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R1_host_removed.fastq.gz"
    R2_HOST_REMOVED="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R2_host_removed.fastq.gz"
    HOST_STATS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_host_removal_stats.txt"

    # Kraken2 output files
    KRAKEN2_OUTPUT="${KRAKEN2_DIR}/kraken2_output_${SAMPLE_ID}"
    KRAKEN2_CLASSIFICATIONS="${KRAKEN2_OUTPUT}/kraken2_output.txt"
    KRAKEN2_REPORT="${KRAKEN2_OUTPUT}/kraken2_report.txt"

    # DIAMOND输出文件
    DIAMOND_OUTPUT="${DIAMOND_DIR}/${SAMPLE_ID}_combined_methane_hits.txt"
    DIAMOND_SUMMARY="${DIAMOND_DIR}/${SAMPLE_ID}_methane_summary.txt"
    
    # Log file
    :

################################################################################
# STEP 1: fastp Quality Control
################################################################################
echo "=========================================="
echo "STEP 1: Quality Control with fastp"
echo "Sample: ${SAMPLE_ID}"
echo "=========================================="

if [ ! -f "${R1_RAW}" ] || [ -z "${R1_RAW}" ] || [ ! -f "${R2_RAW}" ] || [ -z "${R2_RAW}" ]; then
    echo "ERROR: Raw data files not found for sample ${SAMPLE_ID}"
    echo "  Expected: ${RAW_DATA_DIR}/${SAMPLE_ID}_R1*.fastq.gz and R2*.fastq.gz"
    return 1
fi

echo "Found input files:"
echo "  R1: $(basename ${R1_RAW})"
echo "  R2: $(basename ${R2_RAW})"
echo ""

echo "Running fastp..."

# Run fastp for quality control and adapter trimming
# fastp automatically detects and removes adapters
fastp \
    -i ${R1_RAW} \
    -I ${R2_RAW} \
    -o ${R1_CLEAN} \
    -O ${R2_CLEAN} \
    --dedup \
    --dup_calc_accuracy 3 \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread ${THREADS} \
    --json ${FASTP_JSON} \
    --html ${FASTP_HTML} \
    2>&1 | tee "${LOG_DIR}/${SAMPLE_ID}_fastp.log"

if [ $? -ne 0 ]; then
    echo "ERROR: fastp failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "fastp completed!"
echo "Clean reads saved to: ${FASTP_DIR}"
echo "QC report: ${FASTP_HTML}"
echo ""

################################################################################
# STEP 1.5: Host DNA Removal with Bowtie2
################################################################################
echo "=========================================="
echo "STEP 1.5: Host DNA Removal with Bowtie2"
echo "=========================================="

# Update the index prefix based on your screenshot
BT2_INDEX="${BASE_DIR}/Data/reference_dbs/host_genomes/plant_host"

if [ -f "${BT2_INDEX}.1.bt2" ]; then
    echo "Bowtie2 index found: ${BT2_INDEX}"
    echo "Removing plant host reads..."

    # Bowtie2 command for host removal
    # --un-conc-gz: saves paired-end reads that DO NOT map to host (the ones we want)
    # -S /dev/null: we don't need the alignment SAM file, saving disk space
    bowtie2 -x "$BT2_INDEX" \
        -1 ${R1_CLEAN} \
        -2 ${R2_CLEAN} \
        --un-conc-gz "${HOST_REMOVED_DIR}/${SAMPLE_ID}_R%_host_removed.fastq.gz" \
        --threads ${THREADS} \
        --very-sensitive-local \
        -S /dev/null

    # Update paths for downstream analysis
    R1_FOR_ANALYSIS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R1_host_removed.fastq.gz"
    R2_FOR_ANALYSIS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R2_host_removed.fastq.gz"
    
    echo "Host removal completed! Cleaned reads saved to: ${HOST_REMOVED_DIR}"
else
    echo "WARNING: Bowtie2 index NOT found at ${BT2_INDEX}. Skipping host removal."
    R1_FOR_ANALYSIS=${R1_CLEAN}
    R2_FOR_ANALYSIS=${R2_CLEAN}
fi

################################################################################
# STEP 2: Kraken2 Taxonomic Profiling
################################################################################
echo "=========================================="
echo "STEP 2: Taxonomic Profiling with Kraken2"
echo "=========================================="

# Locate kraken2 - check current environment first, then fallback to kraken_env
KRAKEN2_BIN=""
if command -v kraken2 &> /dev/null; then
    KRAKEN2_BIN="kraken2"
    echo "✓ Using kraken2 from current environment: $(which kraken2)"
else
    # Fallback: try to find kraken2 in separate kraken_env (legacy setup)
    KRAKEN2_ENV_PATH="${CONDA_BASE}/envs/kraken_env/bin/kraken2"
    if [ -f "$KRAKEN2_ENV_PATH" ]; then
        KRAKEN2_BIN="$KRAKEN2_ENV_PATH"
        echo "⚠ Using kraken2 from separate environment: $KRAKEN2_BIN"
        echo "  (Consider installing kraken2 in quick_search for better performance)"
    else
        echo "ERROR: kraken2 not found in current environment or kraken_env"
        echo "  Install kraken2 in current environment (recommended):"
        echo "    conda activate quick_search"
        echo "    conda install -c bioconda kraken2"
        return 1
    fi
fi

# Verify kraken2 works
if ! $KRAKEN2_BIN --version &> /dev/null; then
    echo "ERROR: kraken2 binary found but not executable"
    return 1
fi
echo "  Version: $($KRAKEN2_BIN --version 2>&1 | head -1)"

# Set environment variables
export KRAKEN2_DB_PATH=${KRAKEN2_DB}

# Set temporary directory for Kraken2 (avoid /tmp issues on HPC)
if [ -z "$TMPDIR" ]; then
    TMPDIR="${BASE_DIR}/Data/temp"
    mkdir -p "${TMPDIR}"
    export TMPDIR
    export TEMP="${TMPDIR}"
    export TMP="${TMPDIR}"
fi

# Create a unique temporary directory for this sample
SAMPLE_TMPDIR="${TMPDIR}/kraken2_${SAMPLE_ID}_$$"
mkdir -p "${SAMPLE_TMPDIR}"
export TMPDIR="${SAMPLE_TMPDIR}"

# Set Python environment to keep temp files longer
export PYTHONDONTWRITEBYTECODE=1

echo "Using temporary directory: ${TMPDIR}"
# Use more threads for Kraken2/DIAMOND to improve performance
# Only cap on low-memory systems
KRAKEN2_THREADS=${KRAKEN2_THREADS:-${THREADS}}
export OMP_NUM_THREADS=1
# Try to increase stack size but ignore permission warnings
ulimit -s 65536 2>/dev/null || true

# Create output directory with error checking
echo "Creating Kraken2 output directory: ${KRAKEN2_OUTPUT}"
if ! mkdir -p "${KRAKEN2_OUTPUT}"; then
    echo "ERROR: Failed to create Kraken2 output directory: ${KRAKEN2_OUTPUT}"
    echo "Check permissions and disk space."
    return 1
fi

# Verify directory was created
if [ ! -d "${KRAKEN2_OUTPUT}" ]; then
    echo "ERROR: Kraken2 output directory does not exist after creation: ${KRAKEN2_OUTPUT}"
    return 1
fi

echo "Running Kraken2..."
$KRAKEN2_BIN --db ${KRAKEN2_DB} \
    --paired ${R1_FOR_ANALYSIS} ${R2_FOR_ANALYSIS} \
    --output ${KRAKEN2_OUTPUT}/kraken2_output.txt \
    --report ${KRAKEN2_OUTPUT}/kraken2_report.txt \
    --threads ${KRAKEN2_THREADS} \
    --use-names

if [ $? -ne 0 ]; then
    echo "ERROR: Kraken2 classification failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "Kraken2 completed!"
echo "Classification output: ${KRAKEN2_OUTPUT}/kraken2_output.txt"
echo "Taxonomic report: ${KRAKEN2_OUTPUT}/kraken2_report.txt"

# Clean up temporary directory for this sample
rm -rf "${SAMPLE_TMPDIR}" 2>/dev/null || true

echo ""

################################################################################
# STEP 3: DIAMOND Functional Gene Search
################################################################################
echo "=========================================="
echo "STEP 3: Methane Gene Search with DIAMOND"
echo "=========================================="

echo "Using diamond from: $(which diamond)"
echo "Searching for methane metabolism genes..."
# Use a dedicated writable temp dir for DIAMOND
DIAMOND_TMPDIR="${TMPDIR}/diamond_${SAMPLE_ID}_$$"
mkdir -p "${DIAMOND_TMPDIR}" 2>/dev/null || true

gunzip -c ${R1_FOR_ANALYSIS} ${R2_FOR_ANALYSIS} | \
    DIAMOND_TMPDIR="${DIAMOND_TMPDIR}" diamond blastx \
    -d ${METHANE_DB} \
    -q - \
    -o ${DIAMOND_OUTPUT} \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads ${THREADS} \
    --sensitive \
    --query-cover 80 \
    --min-score 40 \
    --max-target-seqs 5 \
    --evalue 1e-5

# Clean DIAMOND temp dir
rm -rf "${DIAMOND_TMPDIR}" 2>/dev/null || true

if [ $? -ne 0 ]; then
    echo "ERROR: DIAMOND search failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "DIAMOND search completed!"
echo "Results saved to: ${DIAMOND_OUTPUT}"
echo ""

################################################################################
# STEP 4: Generate Summary Report
################################################################################
echo "=========================================="
echo "STEP 4: Generating Summary Report"
echo "=========================================="

# Count results
TOTAL_HITS=$(wc -l < ${DIAMOND_OUTPUT})
MMO_HITS=$(grep -ic "methane monooxygenase" ${DIAMOND_OUTPUT})
PMO_HITS=$(grep -iEc "PmoA|pMMO" ${DIAMOND_OUTPUT})
TOTAL_METHANE=$((MMO_HITS + PMO_HITS))

# Generate report
cat > ${DIAMOND_SUMMARY} << EOF
===================================================================
METHANE METABOLISM GENE SEARCH RESULTS - Sample ${SAMPLE_ID}
===================================================================
Analysis Date: $(date +%Y-%m-%d)
Search Tool: DIAMOND blastx
Database: methane_master_db
Search Parameters:
  - Query coverage: >80%
  - Bit score: >40
  - E-value: <1e-5
  - Max target sequences: 5
  - Mode: sensitive

===================================================================
RESULTS SUMMARY
===================================================================

Combined Reads (R1 + R2 from fastp_cleaned):
  Total alignments: ${TOTAL_HITS}
  
  Methane-specific genes:
    - Methane monooxygenase (MMO): ${MMO_HITS} hits
    - PmoA family proteins (pMMO): ${PMO_HITS} hits
    - Total methanotroph genes: ${TOTAL_METHANE} hits

===================================================================
BIOLOGICAL INTERPRETATION
===================================================================

EOF

if [ ${TOTAL_METHANE} -gt 0 ]; then
    cat >> ${DIAMOND_SUMMARY} << EOF
METHANOTROPHS DETECTED
  Sample ${SAMPLE_ID} contains aerobic methane-oxidizing bacteria.
  Key enzymes detected:
  - pMMO (particulate methane monooxygenase): Primary methane oxidation
  - MMO (soluble methane monooxygenase): Alternative pathway
  
  These organisms consume methane (CH4) as carbon and energy source.

===================================================================
ECOLOGICAL SIGNIFICANCE
===================================================================

The presence of methanotrophs suggests:
1. Oxic environment (aerobic methane oxidation)
2. Potential methane sink rather than source
3. Common in wetlands, rice paddies, landfill covers

===================================================================
EOF
else
    cat >> ${DIAMOND_SUMMARY} << EOF
NO SIGNIFICANT METHANE METABOLISM GENES DETECTED
  
  This sample shows low or no methanotroph/methanogen activity.

===================================================================
EOF
fi

echo "Summary report saved to: ${DIAMOND_SUMMARY}"
echo ""

################################################################################
# STEP 5: Gene Abundance Quantification and Normalization
################################################################################
echo "=========================================="
echo "STEP 5: Gene Abundance Quantification"
echo "=========================================="

# Output files for quantification
GENE_COUNTS="${DIAMOND_DIR}/${SAMPLE_ID}_gene_counts.txt"
GENE_RPKM="${DIAMOND_DIR}/${SAMPLE_ID}_gene_rpkm.txt"
FILTERED_HITS="${DIAMOND_DIR}/${SAMPLE_ID}_filtered_methane_hits.txt"

# Filter high-quality hits only (remove ribosomal proteins and low-quality hits)
echo "Filtering high-quality methane gene hits..."
grep -iE "methane|pmo|mmo" ${DIAMOND_OUTPUT} | \
    grep -viE "ribosom|ribokinase|ribonucleoside|ammonia" > ${FILTERED_HITS}

FILTERED_COUNT=$(wc -l < ${FILTERED_HITS})
echo "Filtered hits: ${FILTERED_COUNT} (removed non-methane genes)"

# Count reads per gene
echo "Counting reads per gene..."
awk '{print $7}' ${FILTERED_HITS} | sort | uniq -c | sort -rn > ${GENE_COUNTS}

echo "Top 10 most abundant genes:"
head -10 ${GENE_COUNTS}
echo ""

# Calculate normalization factors
# Get total mapped reads (from original DIAMOND output)
TOTAL_MAPPED_READS=$(awk '{print $1}' ${DIAMOND_OUTPUT} | sort -u | wc -l)

# Get total sequencing reads (from both R1 and R2) - optimized single pass
TOTAL_READS=$(( $(zcat ${R1_FOR_ANALYSIS} ${R2_FOR_ANALYSIS} | wc -l) / 4 ))
# Calculate millions using bash arithmetic (avoids bc dependency)
TOTAL_READS_MILLIONS=$(awk "BEGIN {printf \"%.4f\", ${TOTAL_READS} / 1000000}")

echo "Normalization factors:"
echo "  Total reads (R1+R2): ${TOTAL_READS}"
echo "  Total reads in millions: ${TOTAL_READS_MILLIONS}M"
echo "  Total mapped reads: ${TOTAL_MAPPED_READS}"
echo ""

# Calculate RPKM for each gene
# RPKM = (Read counts * 10^9) / (Gene length in bp * Total mapped reads)
# For protein-level alignment, we use alignment length as proxy for gene length
echo "Calculating RPKM values..."

awk -v total_reads="${TOTAL_READS_MILLIONS}" -v filtered_hits="${FILTERED_HITS}" '
BEGIN {
    print "Gene_Name\tRead_Count\tAvg_Length\tRPKM"
}
{
    gene_name = $0
    gsub(/^[[:space:]]+/, "", gene_name)
    split(gene_name, arr, " ")
    count = arr[1]
    gene = ""
    for(i=2; i<=length(arr); i++) gene = gene " " arr[i]
    gsub(/^[[:space:]]+/, "", gene)
    
    # Get average alignment length for this gene
    cmd = "grep -F \"" gene "\" " filtered_hits " | awk '"'"'{sum+=$4; n++} END {if(n>0) print sum/n; else print 1000}'"'"'"
    cmd | getline avg_length
    close(cmd)
    
    # Calculate RPKM: (count * 1000) / (avg_length * total_reads_millions)
    rpkm = (count * 1000) / (avg_length * total_reads)
    
    printf "%s\t%d\t%.1f\t%.4f\n", gene, count, avg_length, rpkm
}
' ${GENE_COUNTS} > ${GENE_RPKM}

echo "RPKM calculation completed!"
echo "Results saved to: ${GENE_RPKM}"
echo ""

echo "Top 10 genes by RPKM:"
head -11 ${GENE_RPKM} | tail -10
echo ""

################################################################################
# STEP 6: Display Results
################################################################################
echo "=========================================="
echo "ANALYSIS COMPLETE - Sample ${SAMPLE_ID}"
echo "=========================================="
echo ""
echo "Output Files:"
echo "  fastp cleaned reads:"
echo "    - ${R1_CLEAN}"
echo "    - ${R2_CLEAN}"
echo "    - QC report: ${FASTP_HTML}"
if [[ -f "${HOST_REF}" ]]; then
    echo "  Host-removed reads:"
    echo "    - ${R1_HOST_REMOVED}"
    echo "    - ${R2_HOST_REMOVED}"
    echo "    - Host removal stats: ${HOST_STATS}"
fi
echo ""
echo "  Kraken2 results:"
echo "    - ${KRAKEN2_OUTPUT}/kraken2_output.txt"
echo "    - ${KRAKEN2_OUTPUT}/kraken2_report.txt"
echo ""
echo "  DIAMOND results:"
echo "    - ${DIAMOND_OUTPUT}"
echo "    - ${DIAMOND_SUMMARY}"
echo ""
echo "  Gene abundance analysis:"
echo "    - ${FILTERED_HITS}"
echo "    - ${GENE_COUNTS}"
echo "    - ${GENE_RPKM}"
echo ""
echo "Quick Results:"
echo "  Total DIAMOND hits: ${TOTAL_HITS}"
echo "  Filtered methane gene hits: ${FILTERED_COUNT}"
echo "  Methane monooxygenase: ${MMO_HITS} hits"
echo "  PmoA family proteins: ${PMO_HITS} hits"
echo "  Total methanotroph genes: ${TOTAL_METHANE} hits"
echo "  Total sequencing reads: ${TOTAL_READS}"
echo ""
    local END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))
    echo "All steps completed successfully!"
    echo "Time elapsed: $((ELAPSED / 60)) minutes $((ELAPSED % 60)) seconds"
    echo "=========================================="
    echo ""
    ) > >(tee -a "$LOG_FILE") 2>&1
    return $?
}

################################################################################
# Main execution loop
################################################################################
MAIN_START_TIME=$(date +%s)
SUCCESS_COUNT=0
FAIL_COUNT=0
FAILED_SAMPLES=()

# Parallel processing function
if [ ${MAX_PARALLEL_SAMPLES} -gt 1 ]; then
    echo "Processing ${#SAMPLE_LIST[@]} samples with ${MAX_PARALLEL_SAMPLES} parallel jobs..."
    echo ""
    
    # Use GNU parallel or xargs for parallel processing
    if command -v parallel &> /dev/null; then
        export -f process_sample
        export BASE_DIR THREADS RAW_DATA_DIR FASTP_DIR KRAKEN2_DIR DIAMOND_DIR LOG_DIR
        export KRAKEN2_DB METHANE_DB KRAKEN2_THREADS RAW_DATA_SUBDIR
        
        printf '%s\n' "${SAMPLE_LIST[@]}" | \
            parallel -j ${MAX_PARALLEL_SAMPLES} --line-buffer \
            'if process_sample {}; then echo "✓ {} completed"; else echo "✗ {} failed"; fi'
        
        # Count successes/failures from logs
        for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
            if grep -q "All steps completed successfully" "${LOG_DIR}/${SAMPLE_ID}_quick_search.log" 2>/dev/null; then
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                FAIL_COUNT=$((FAIL_COUNT + 1))
                FAILED_SAMPLES+=("$SAMPLE_ID")
            fi
        done
    else
        echo "WARNING: GNU parallel not found, falling back to sequential processing"
        echo "Install with: conda install -c conda-forge parallel"
        echo ""
        MAX_PARALLEL_SAMPLES=1
    fi
fi

# Sequential processing (default or fallback)
if [ ${MAX_PARALLEL_SAMPLES} -eq 1 ]; then
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo ""
        echo "################################################################################"
        echo "# Processing sample ${SAMPLE_ID} ($(($SUCCESS_COUNT + $FAIL_COUNT + 1))/${#SAMPLE_LIST[@]})"
        echo "################################################################################"
        
        if process_sample "$SAMPLE_ID"; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            echo "✓ Sample ${SAMPLE_ID} completed successfully"
        else
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAILED_SAMPLES+=("$SAMPLE_ID")
            echo "✗ Sample ${SAMPLE_ID} failed"
        fi
    done
fi

################################################################################
# Final Summary
################################################################################
MAIN_END_TIME=$(date +%s)
TOTAL_ELAPSED=$((MAIN_END_TIME - MAIN_START_TIME))

echo ""
echo "========================================================================"
echo "PIPELINE SUMMARY"
echo "========================================================================"
echo "Total samples processed: ${#SAMPLE_LIST[@]}"
echo "Successful: ${SUCCESS_COUNT}"
echo "Failed: ${FAIL_COUNT}"
echo "Total time: $((TOTAL_ELAPSED / 3600))h $((TOTAL_ELAPSED % 3600 / 60))m $((TOTAL_ELAPSED % 60))s"
echo ""

if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "Failed samples:"
    for SAMPLE in "${FAILED_SAMPLES[@]}"; do
        echo "  - ${SAMPLE}"
    done
    echo ""
    echo "Check logs in: ${LOG_DIR}/"
    echo "========================================================================"
    exit 1
fi

echo "All samples completed successfully!"
echo ""
echo "Output locations:"
echo "  Quality-controlled reads: Data/processed_data/fastp_cleaned/"
echo "  Taxonomic profiles:       Data/processed_data/kraken2_output/"
echo "  Gene detection results:   Data/functional_analysis/methane_genes/"
echo "  Logs:                     Logs/"
echo "======================================================================="
