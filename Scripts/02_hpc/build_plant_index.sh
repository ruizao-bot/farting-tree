#!/bin/bash
#PBS -N build_plant_index
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=02:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
source /rds/general/user/jc224/home/miniconda3/etc/profile.d/conda.sh
conda activate kraken_env

bowtie2-build --threads 8 combined_host.fasta plant_host