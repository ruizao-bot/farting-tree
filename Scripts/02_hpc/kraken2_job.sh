#!/bin/bash
#PBS -N kraken2_job
#PBS -l select=1:ncpus=24:mem=150gb
#PBS -l walltime=24:00:00
#PBS -j oe

# 进入脚本所在的目录
cd $PBS_O_WORKDIR


source /rds/general/user/jc224/home/miniconda3/etc/profile.d/conda.sh
conda activate kraken_env


DB_DIR="/rds/general/user/jc224/home/meta/Data/reference_dbs"
DATA_DIR="/rds/general/user/jc224/home/meta/Data/processed_data/host_removed"
OUTDIR="/rds/general/user/jc224/home/meta/Results/kraken_host_removal"
mkdir -p "$OUTDIR"

for fq1 in ${DATA_DIR}/*_R1*.fastq.gz; do
    sample=$(basename "$fq1" | sed 's/_R1.*//')
    fq2="${DATA_DIR}/${sample}_R2.fastq.gz"
    
    # Handle alternative R2 naming patterns
    if [ ! -f "$fq2" ]; then
        fq2=$(echo "$fq1" | sed 's/_R1/_R2/')
    fi
    
    if [ -f "$fq2" ]; then
        echo "Processing $sample (host-removed data)..."
        kraken2 --db $DB_DIR \
            --threads 24 \
            --report "$OUTDIR/${sample}.report" \
            --output "$OUTDIR/${sample}.kraken" \
            --paired --gzip-compressed "$fq1" "$fq2"
    else
        echo "WARNING: R2 file not found for $sample, skipping..."
    fi
done
