#!/bin/bash -l
#SBATCH --job-name=shotgun_analysis  # Job name
#SBATCH --output=shotgun_analysis_%j.out  # Standard output log
#SBATCH --error=shotgun_analysis_%j.err   # Standard error log
#SBATCH --ntasks=1                       # Number of tasks (processes)
#SBATCH --cpus-per-task=16               # Number of CPU cores per task
#SBATCH --mem=64G                        # Memory per node
#SBATCH --time=48:00:00                  # Time limit (48 hours)
#SBATCH --partition=large_336              # Partition name

# Initialize Conda
eval "$(conda shell.bash hook)"
conda activate qiime2-moshpit

# Navigate to the project directory
cd /home/jiayi/QIIME/Scripts

# Run the shotgun analysis script
bash shotgun.sh

# Deactivate Conda environment
conda deactivate