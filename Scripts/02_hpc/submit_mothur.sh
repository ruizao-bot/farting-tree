#!/bin/bash
#PBS -N mothur_pmoA
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd "$PBS_O_WORKDIR"

# Make sure conda is available
source ~/.bashrc
conda activate bio_env

# Run mothur with your batch file
mothur ../01_pipelines/run_pmoA.mth