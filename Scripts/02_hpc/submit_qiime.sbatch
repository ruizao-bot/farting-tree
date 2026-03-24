#!/bin/bash
#SBATCH --job-name=qiime_process          # Job name
#SBATCH --output=Logs/slurm-%j.out        # Stdout log (%j = Job ID)
#SBATCH --error=Logs/slurm-%j.err         # Stderr log
#SBATCH --time=24:00:00                   # Time limit
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # CPUs per task (kept in sync with N_THREADS)
#SBATCH --mem=64G                         # Memory
# #SBATCH --partition=standard            # Partition/queue (comment out or check with 'sinfo')

# Minimal Slurm submission script for QIIME pipeline.
# Submit from the repo root on the HPC node:
#   sbatch Scripts/submit_qiime.sbatch

set -euo pipefail

# Move to repository root to ensure relative paths resolve correctly
cd "$(dirname "$0")/.."

# Load conda (choose what matches your cluster)
# Option A: environment modules (typical on academic clusters)
if command -v module >/dev/null 2>&1; then
  module load anaconda3 2>/dev/null || module load miniconda3 2>/dev/null || true
fi
# Option B: source from a known install if modules are unavailable
if ! command -v conda >/dev/null 2>&1; then
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  fi
else
  # Ensure conda activation script is sourced
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi


# Activate QIIME environment (change name if needed)
conda activate qiime

# Ensure Logs directory exists for Slurm outputs
mkdir -p Logs

echo "Job started on $(hostname) at $(date)"

# Pass allocated CPUs to the pipeline
export N_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Interactive run (prompts for input):
# bash Scripts/main.sh

# Non-interactive run (recommended on HPC, no prompts):
bash Scripts/main.sh --non-interactive

echo "Job finished at $(date)"
