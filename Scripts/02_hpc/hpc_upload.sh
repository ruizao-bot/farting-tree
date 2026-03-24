#!/bin/bash

# Upload QIIME pipeline to HPC cluster and submit as Slurm job
# Usage: bash hpc_upload.sh HPC_HOST HPC_PROJECT_DIR [SLURM_OPTIONS]
# Example: bash hpc_upload.sh user@hpc.edu /scratch/user/QIIME --time=48:00:00 --mem=128G

set -euo pipefail

# ============================================================================
# CUSTOMIZE THESE SETTINGS FOR YOUR HPC ENVIRONMENT
# ============================================================================
CONDA_ENV_NAME="qiime"           # Name of your QIIME2 conda environment on HPC
CONDA_INIT_METHOD="auto"         # "auto", "conda_base", or "module"
                                 # auto: tries to find conda automatically
                                 # conda_base: sources conda from $(conda info --base)
                                 # module: uses 'module load anaconda' or 'module load miniconda'

# Slurm default parameters (override with command-line options)
SLURM_TIME="24:00:00"
SLURM_CPUS="8"
SLURM_MEM="64G"
SLURM_PARTITION="standard"
SLURM_JOB_NAME="qiime_pipeline"

# Configuration
HPC_HOST="jiayi@harvey.dept.imperial.ac.uk"           # Your HPC login
HPC_PROJECT_DIR="/home/jiayi/QIIME"                    # Your project path on HPC
LOCAL_PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# Help message
show_help() {
    cat << 'EOF'
Usage: bash hpc_upload.sh HPC_HOST HPC_PROJECT_DIR [OPTIONS]

Upload QIIME pipeline to HPC cluster and submit as Slurm job.

Required Arguments:
  HPC_HOST          SSH connection string (e.g., user@hpc.university.edu)
  HPC_PROJECT_DIR   Remote project directory on HPC (e.g., /scratch/user/QIIME)

Slurm Options (all optional):
  --time=TIME       Job time limit (default: 24:00:00)
  --cpus=N          Number of CPUs (default: 8)
  --mem=SIZE        Memory allocation (default: 64G)
  --partition=NAME  Slurm partition (default: standard)
  --job-name=NAME   Job name (default: qiime_pipeline)

Examples:
  # Basic usage with defaults
  bash hpc_upload.sh user@hpc.edu /scratch/user/QIIME

  # Custom resources for large dataset
  bash hpc_upload.sh user@hpc.edu /scratch/user/QIIME \
    --time=48:00:00 --mem=128G --cpus=16

  # Different partition
  bash hpc_upload.sh user@hpc.edu /scratch/user/QIIME \
    --partition=large --time=72:00:00

Steps performed:
  1. Syncs local project directory to HPC using rsync
  2. Creates sbatch script with specified Slurm parameters
  3. Submits the job to Slurm queue
  4. Shows job status and monitoring commands

Before running:
  - Set up SSH key authentication: ssh-copy-id user@hpc.edu
  - Customize CONDA_ENV_NAME and CONDA_INIT_METHOD at top of this script
  - Ensure QIIME2 conda environment exists on HPC

Pipeline Parameters:
  To run in non-interactive mode, set environment variables before uploading:
    export TRUNC_LEN_F=250
    export TRUNC_LEN_R=250
    export N_THREADS=16
    export DECONTAMINATION_CHOICE=1
    export DECONTAM_THRESHOLD=0.1
    export PRIMER_CHOICE=3
    export CLASSIFIER_CHOICE=1
  
  Or run interactively (main.sh will prompt for parameters)

EOF
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_help
            ;;
        --time=*)
            SLURM_TIME="${1#*=}"
            shift
            ;;
        --cpus=*)
            SLURM_CPUS="${1#*=}"
            shift
            ;;
        --mem=*)
            SLURM_MEM="${1#*=}"
            shift
            ;;
        --partition=*)
            SLURM_PARTITION="${1#*=}"
            shift
            ;;
        --job-name=*)
            SLURM_JOB_NAME="${1#*=}"
            shift
            ;;
        *)
            if [[ -z "${HPC_HOST}" ]]; then
                HPC_HOST="$1"
            elif [[ -z "${HPC_PROJECT_DIR}" ]]; then
                HPC_PROJECT_DIR="$1"
            else
                echo "ERROR: Unknown argument: $1"
                echo "Run with --help for usage information"
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate required arguments
if [[ -z "${HPC_HOST}" ]] || [[ -z "${HPC_PROJECT_DIR}" ]]; then
    echo "ERROR: Both HPC_HOST and HPC_PROJECT_DIR are required"
    echo "Run with --help for usage information"
    exit 1
fi

echo "=== QIIME Pipeline HPC Upload & Submit ==="
echo "Local project:  ${LOCAL_PROJECT_DIR}"
echo "HPC host:       ${HPC_HOST}"
echo "HPC directory:  ${HPC_PROJECT_DIR}"
echo ""
echo "Slurm Configuration:"
echo "  Job name:     ${SLURM_JOB_NAME}"
echo "  Time limit:   ${SLURM_TIME}"
echo "  CPUs:         ${SLURM_CPUS}"
echo "  Memory:       ${SLURM_MEM}"
echo "  Partition:    ${SLURM_PARTITION}"
echo ""

# Test SSH connection
echo "Testing SSH connection..."
if ! ssh -o ConnectTimeout=10 "${HPC_HOST}" "echo 'Connection successful'" 2>/dev/null; then
    echo "ERROR: Cannot connect to ${HPC_HOST}"
    echo "Please check:"
    echo "  1. SSH key is set up (ssh-copy-id ${HPC_HOST})"
    echo "  2. Host is correct and accessible"
    exit 1
fi
echo "✓ SSH connection successful"
echo ""

# Create remote directory if needed
echo "Creating remote directory..."
ssh "${HPC_HOST}" "mkdir -p ${HPC_PROJECT_DIR}"
echo "✓ Remote directory ready"
echo ""

# Sync project files to HPC
echo "Syncing project files to HPC..."
rsync -avz --progress \
    --exclude '.git' \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    --exclude '.DS_Store' \
    --exclude 'Results/' \
    --exclude 'Data/processed_data/' \
    --exclude 'Data/temp/' \
    --exclude 'Logs/' \
    "${LOCAL_PROJECT_DIR}/" \
    "${HPC_HOST}:${HPC_PROJECT_DIR}/"
echo "✓ Files synced"
echo ""

# Create sbatch script on HPC
echo "Creating sbatch job script on HPC..."
ssh "${HPC_HOST}" "cat > ${HPC_PROJECT_DIR}/submit_pipeline.sbatch" << EOF
#!/bin/bash
#SBATCH --job-name=${SLURM_JOB_NAME}
#SBATCH --output=${HPC_PROJECT_DIR}/Logs/slurm-pipeline-%j.out
#SBATCH --error=${HPC_PROJECT_DIR}/Logs/slurm-pipeline-%j.err
#SBATCH --time=${SLURM_TIME}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --partition=${SLURM_PARTITION}

# Auto-generated sbatch script for QIIME pipeline
# Generated by hpc_upload.sh on $(date)

set -euo pipefail

# Change to project directory
cd ${HPC_PROJECT_DIR}

# Create Logs directory
mkdir -p Logs

# Initialize conda based on CONDA_INIT_METHOD setting
case "${CONDA_INIT_METHOD}" in
    conda_base)
        # Source conda from base installation
        CONDA_BASE="\$(conda info --base 2>/dev/null || echo "\${HOME}/miniconda3")"
        if [[ -f "\${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
            source "\${CONDA_BASE}/etc/profile.d/conda.sh"
        else
            echo "ERROR: Cannot find conda.sh at \${CONDA_BASE}/etc/profile.d/conda.sh"
            exit 1
        fi
        ;;
    module)
        # Use environment modules
        if command -v module >/dev/null 2>&1; then
            module load anaconda 2>/dev/null || module load miniconda 2>/dev/null || module load conda 2>/dev/null || true
        else
            echo "ERROR: module command not found"
            exit 1
        fi
        ;;
    auto|*)
        # Try to auto-detect conda
        if [[ -f "\$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]]; then
            source "\$(conda info --base)/etc/profile.d/conda.sh"
        elif command -v module >/dev/null 2>&1; then
            module load anaconda 2>/dev/null || module load miniconda 2>/dev/null || true
        elif [[ -f "\${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
            source "\${HOME}/miniconda3/etc/profile.d/conda.sh"
        elif [[ -f "\${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
            source "\${HOME}/anaconda3/etc/profile.d/conda.sh"
        else
            echo "ERROR: Cannot find conda. Please set CONDA_INIT_METHOD in hpc_upload.sh"
            exit 1
        fi
        ;;
esac

# Activate QIIME environment
echo "Activating conda environment: ${CONDA_ENV_NAME}"
conda activate ${CONDA_ENV_NAME}

# Verify QIIME2 is available
if ! command -v qiime >/dev/null 2>&1; then
    echo "ERROR: qiime command not found after activating ${CONDA_ENV_NAME}"
    exit 1
fi

echo "QIIME2 version:"
qiime --version

# Set number of threads to match allocated CPUs
export N_THREADS=${SLURM_CPUS}

# Run pipeline in non-interactive mode
# Note: Set other environment variables (TRUNC_LEN_F, etc.) before running hpc_upload.sh
# or the pipeline will use defaults
export NON_INTERACTIVE=true
bash Scripts/01_pipelines/16S.sh --non-interactive

echo "Pipeline completed at \$(date)"
EOF

echo "✓ sbatch script created"
echo ""

# Submit job to Slurm
echo "Submitting job to Slurm queue..."
JOB_ID=$(ssh "${HPC_HOST}" "cd ${HPC_PROJECT_DIR} && sbatch submit_pipeline.sbatch" | grep -oP '\d+')

if [[ -n "${JOB_ID}" ]]; then
    echo "✓ Job submitted successfully!"
    echo "  Job ID: ${JOB_ID}"
    echo ""
    echo "Monitor job status with:"
    echo "  ssh ${HPC_HOST} 'squeue -j ${JOB_ID}'"
    echo ""
    echo "View job output:"
    echo "  ssh ${HPC_HOST} 'tail -f ${HPC_PROJECT_DIR}/Logs/slurm-pipeline-${JOB_ID}.out'"
    echo ""
    echo "Cancel job if needed:"
    echo "  ssh ${HPC_HOST} 'scancel ${JOB_ID}'"
else
    echo "ERROR: Job submission failed"
    exit 1
fi
