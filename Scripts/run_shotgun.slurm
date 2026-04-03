#!/bin/bash
#SBATCH -J shotgun_quick       # Job name
#SBATCH -p large_336           # Partition (adjust based on your cluster)
#SBATCH -N 1                   # Number of nodes
#SBATCH --ntasks-per-node=16   # CPUs per node
#SBATCH --mem=128G             # Memory (Kraken2 needs substantial RAM)
#SBATCH -t 48:00:00            # Time limit (48 hours)
#SBATCH -o Logs/slurm_%j.out   # Standard output log
#SBATCH -e Logs/slurm_%j.err   # Standard error log

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Record start time
echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "=========================================="
echo ""

# Navigate to project directory
cd /home/jc224/shotgun

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate quick_search

# Environment check
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo ""

# ============================================================================
# Configuration - Adjust these parameters as needed
# ============================================================================

# Number of threads for parallel tools (fastp, DIAMOND, Kraken2, Bowtie2)
export THREADS=${THREADS:-16}

# Maximum number of samples to process in parallel
# Set to 1 for sequential processing (safer for large samples)
# Set higher if samples are small and you have GNU parallel installed
export MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-1}

# Force re-run all steps even if outputs exist
# export FORCE=true

# Raw data subdirectory name (default: shotgun)
# export RAW_DATA_SUBDIR=shotgun

# ============================================================================
# Sample selection - Choose ONE of the following methods
# ============================================================================

# METHOD 1: Auto-detect all samples in Data/raw_data/shotgun/
AUTO_MODE="--auto"

# METHOD 2: Specify individual sample IDs
# SAMPLES="53394 53395 53396"
# AUTO_MODE=""

# METHOD 3: Use a sample list file
# Create a file with one sample ID per line, e.g., Config/samples.txt
# SAMPLE_FILE="Config/samples.txt"
# AUTO_MODE="--sample-list $SAMPLE_FILE"

# ============================================================================
# Run the pipeline
# ============================================================================

echo "Starting shotgun metagenomics pipeline..."
echo "Threads: $THREADS"
echo "Max parallel samples: $MAX_PARALLEL_SAMPLES"
echo ""

if [ -n "${AUTO_MODE:-}" ]; then
    # Auto-detect mode
    bash Scripts/shotgun_quick.sh \
        $AUTO_MODE \
        --threads "$THREADS" \
        --parallel "$MAX_PARALLEL_SAMPLES"
else
    # Manual sample list
    bash Scripts/shotgun_quick.sh \
        $SAMPLES \
        --threads "$THREADS" \
        --parallel "$MAX_PARALLEL_SAMPLES"
fi

# ============================================================================
# Job completion
# ============================================================================

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job finished at: $(date)"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Status: SUCCESS ✓"
else
    echo "Status: FAILED (exit code: $EXIT_CODE)"
fi
echo "=========================================="

exit $EXIT_CODE
