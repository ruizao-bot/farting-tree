#!/bin/bash -l
#SBATCH --job-name=quick_search        # Job name
#SBATCH --output=quick_search_%j.out   # Standard output log (%j = job ID)
#SBATCH --error=quick_search_%j.err    # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=64G                      # Memory per node
#SBATCH --time=24:00:00                # Time limit (24 hours)
#SBATCH --partition=large_336          # Partition name

################################################################################
# SLURM Job Script for Quick Metagenomic Search Pipeline
# 
# Usage on HPC:
#   Single sample:
#     sbatch submit_quick_search.sbatch 53394
#   
#   Multiple samples:
#     sbatch submit_quick_search.sbatch 53394 53395 53396
#   
#   From sample list file:
#     sbatch submit_quick_search.sbatch --sample-list samples.txt
#   
#   Auto-detect all samples:
#     sbatch submit_quick_search.sbatch --auto
#
################################################################################

# Print job information
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start Time: $(date)"
echo "=========================================="
echo ""

# Set environment variables for HPC
export BASE_DIR="/home/jiayi/QIIME"
export THREADS=$SLURM_CPUS_PER_TASK
export MAX_PARALLEL_SAMPLES=1  # Set to 2-4 for parallel sample processing
export CONDA_BASE="/home/jiayi/miniconda3"  # Adjust if your conda is elsewhere
export RAW_DATA_SUBDIR="shotgun"  # Specify shotgun subdirectory

# Set temporary directory to avoid /tmp space issues
export TMPDIR="${BASE_DIR}/Data/temp"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

# Initialize Conda
if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
else
    eval "$(conda shell.bash hook)"
fi

# Activate the main environment (bbduk and diamond should be here)
conda activate quick_search

# Set database paths (adjust these for your HPC environment)
export SINGLEM_METAPACKAGE="${BASE_DIR}/Data/reference_dbs/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"
export METHANE_DB="${BASE_DIR}/Data/reference_dbs/DIAMOND/methane_master_db"

# Find BBMap adapters
if [ -f "${CONDA_BASE}/envs/qiime2-moshpit/opt/bbmap*/resources/adapters.fa" ]; then
    export ADAPTERS=$(ls ${CONDA_BASE}/envs/qiime2-moshpit/opt/bbmap*/resources/adapters.fa 2>/dev/null | head -1)
fi

# Navigate to the scripts directory
cd ${BASE_DIR}/Scripts

# Parse arguments: support --sample-list and single sample ID
SAMPLE_LIST_FILE=""
SAMPLE_ARG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample-list)
            SAMPLE_LIST_FILE="$2"
            shift 2
            ;;
        --auto)
            SAMPLE_ARG="--auto"
            shift
            ;;
        *)
            SAMPLE_ARG="$1"
            shift
            ;;
    esac
done

# If a sample list file was provided and we're not running as an array, auto-submit an array job
if [[ -n "$SAMPLE_LIST_FILE" && -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if [ ! -f "$SAMPLE_LIST_FILE" ]; then
        echo "ERROR: Sample list file not found: $SAMPLE_LIST_FILE"
        exit 1
    fi
    NUM_TASKS=$(grep -cvE '^\s*#|^$' "$SAMPLE_LIST_FILE")
    if [ $NUM_TASKS -eq 0 ]; then
        echo "ERROR: Sample list file contains no samples: $SAMPLE_LIST_FILE"
        exit 1
    fi
    echo "No SLURM_ARRAY_TASK_ID detected; submitting array job with $NUM_TASKS tasks..."
    sbatch --array=1-$NUM_TASKS "$0" --sample-list "$SAMPLE_LIST_FILE"
    # Exit current (non-array) submission after spawning array
    exit 0
fi

# Determine which sample to run
if [[ -n "$SAMPLE_LIST_FILE" ]]; then
    # Running as an array task
    SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST_FILE" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [ -z "$SAMPLE_ID" ]; then
        echo "ERROR: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
        exit 1
    fi
elif [[ "$SAMPLE_ARG" == "--auto" ]]; then
    SAMPLE_ID="--auto"
elif [[ -n "$SAMPLE_ARG" ]]; then
    SAMPLE_ID="$SAMPLE_ARG"
else
    echo "ERROR: No sample specified. Provide a sample ID or use --sample-list samples.txt"
    exit 1
fi

# Run the quick search script for a single sample (or in auto mode)
echo "Running quick_search.sh for sample: $SAMPLE_ID"
echo ""

if [[ "$SAMPLE_ID" == "--auto" ]]; then
    bash quick_search.sh --auto
else
    bash quick_search.sh "$SAMPLE_ID"
fi

# Capture exit status
EXIT_STATUS=$?

# Deactivate Conda environment
conda deactivate

echo ""
echo "=========================================="
echo "Job completed with exit status: $EXIT_STATUS"
echo "End Time: $(date)"
echo "=========================================="

exit $EXIT_STATUS
