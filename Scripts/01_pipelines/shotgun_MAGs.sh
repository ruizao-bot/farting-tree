#!/usr/bin/env bash

# QIIME2 Shotgun Data Analysis Pipeline
# This script processes shotgun metagenomic data using QIIME2
# Ensure QIIME2 is installed and the environment is activated before running

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

eval "$(conda shell.bash hook)"
# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/Data"
RAW_DATA_DIR="${DATA_DIR}/raw_data/shotgun"
PROCESSED_DATA_DIR="${DATA_DIR}/processed_data"
REFERENCE_DB_DIR="${DATA_DIR}/reference_dbs"
METADATA_FILE="${DATA_DIR}/metadata/metadata.tsv"
MANIFEST_FILE="${RAW_DATA_DIR}/manifest.tsv"
LOG_FILE="${PROJECT_DIR}/Logs/shotgun_pipeline.log"
ENV_NAME="qiime2-moshpit"
ASSEMBLER="megahit"  # Options: megahit, spades
ASSEMBLY_MODE="per-sample"  # Options: per-sample, coassembly

# Activate qiime2-moshpit environment at the beginning of the script
conda activate qiime2-moshpit

# Create necessary directories
mkdir -p "${PROCESSED_DATA_DIR}"
mkdir -p "${PROJECT_DIR}/Logs"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if output file exists and skip step
check_skip() {
    local output_file="$1"
    local step_name="$2"
    if [[ -f "${output_file}" ]]; then
        log "${step_name}: Output file already exists (${output_file}). Skipping..."
        return 0  # Skip
    fi
    return 1  # Don't skip
}

# This step has been consolidated into step2_quality_control

# Step 1: Import Data
step1_import_data() {
    log "Starting Step 1: Importing shotgun data"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/demux-paired-end.qza" "Step 1"; then
        return 0
    fi

    # Check if manifest file exists
    if [[ ! -f "${MANIFEST_FILE}" ]]; then
        error_exit "Manifest file not found at ${MANIFEST_FILE}. Please create it before running the pipeline."
    fi

    # Import data
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST_FILE}" \
        --output-path "${PROCESSED_DATA_DIR}/demux-paired-end.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to import data."

    log "Data imported successfully. Output: ${PROCESSED_DATA_DIR}/demux-paired-end.qza"
}

# Step 2: Comprehensive Quality Control (FastQC + fastp + MultiQC)
step2_quality_control() {
    log "Starting Step 2: Comprehensive Quality Control"

    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/fastp-cleaned-seqs.qza" "Step 2"; then
        return 0
    fi

    # Set up directories
    QC_DIR="${PROJECT_DIR}/Data/quality_control"
    FASTQC_RAW_DIR="${QC_DIR}/fastqc_raw"
    FASTQC_CLEANED_DIR="${QC_DIR}/fastqc_cleaned"
    FASTP_DIR="${QC_DIR}/fastp_reports"
    MULTIQC_DIR="${QC_DIR}/multiqc_report"
    TEMP_DIR="${PROCESSED_DATA_DIR}/temp_fastp"
    
    mkdir -p "${FASTQC_RAW_DIR}" "${FASTQC_CLEANED_DIR}" "${FASTP_DIR}" "${MULTIQC_DIR}" "${TEMP_DIR}/cleaned"

    # Check manifest file exists
    if [[ ! -f "${MANIFEST_FILE}" ]]; then
        error_exit "Original manifest file not found at ${MANIFEST_FILE}"
    fi

    # ========== 1. FastQC on raw data ==========
    if [[ ! -f "${FASTQC_RAW_DIR}/.done" ]]; then
        log "[QC Step 1/4] Running FastQC on raw data..."
        fastqc -t 4 -o "${FASTQC_RAW_DIR}" "${RAW_DATA_DIR}"/*.fastq.gz || \
            error_exit "Failed to run FastQC on raw data."
        touch "${FASTQC_RAW_DIR}/.done"
        log "FastQC on raw data completed."
    else
        log "[QC Step 1/4] FastQC on raw data already completed. Skipping."
    fi

    # ========== 2. fastp processing ==========
    if [[ ! -f "${FASTP_DIR}/.done" ]]; then
        log "[QC Step 2/4] Running fastp for adapter trimming and quality filtering..."
        
        # Read manifest line by line (skip header)
        while IFS=$'\t' read -r SAMPLE R1 R2; do
            # Skip empty lines
            if [[ -z "${SAMPLE}" ]]; then
                continue
            fi
            
            # Trim whitespace
            SAMPLE=$(echo "${SAMPLE}" | xargs)
            R1=$(echo "${R1}" | xargs)
            R2=$(echo "${R2}" | xargs)
            
            # Skip if already processed
            if [[ -f "${TEMP_DIR}/cleaned/${SAMPLE}_R1.fastq.gz" && -f "${TEMP_DIR}/cleaned/${SAMPLE}_R2.fastq.gz" ]]; then
                log "Sample ${SAMPLE} already processed. Skipping."
                continue
            fi
            
            if [[ ! -f "${R1}" ]]; then
                error_exit "R1 file not found: ${R1} (sample: ${SAMPLE})"
            fi
            if [[ ! -f "${R2}" ]]; then
                error_exit "R2 file not found: ${R2} (sample: ${SAMPLE})"
            fi

            log "Processing sample: ${SAMPLE}"

            fastp \
                -i "${R1}" \
                -I "${R2}" \
                -o "${TEMP_DIR}/cleaned/${SAMPLE}_R1.fastq.gz" \
                -O "${TEMP_DIR}/cleaned/${SAMPLE}_R2.fastq.gz" \
                -h "${FASTP_DIR}/${SAMPLE}_fastp.html" \
                -j "${FASTP_DIR}/${SAMPLE}_fastp.json" \
                -q 20 -u 30 -n 5 -l 50 -w 4 || \
                error_exit "fastp failed for sample ${SAMPLE}"
        done < <(tail -n +2 "${MANIFEST_FILE}")
        touch "${FASTP_DIR}/.done"
        log "fastp processing completed."
    else
        log "[QC Step 2/4] fastp processing already completed. Skipping."
    fi

    # ========== 3. FastQC on cleaned data ==========
    if [[ ! -f "${FASTQC_CLEANED_DIR}/.done" ]]; then
        log "[QC Step 3/4] Running FastQC on cleaned data..."
        fastqc -t 4 -o "${FASTQC_CLEANED_DIR}" "${TEMP_DIR}/cleaned/"*.fastq.gz || \
            error_exit "Failed to run FastQC on cleaned data."
        touch "${FASTQC_CLEANED_DIR}/.done"
        log "FastQC on cleaned data completed."
    else
        log "[QC Step 3/4] FastQC on cleaned data already completed. Skipping."
    fi

    # ========== 4. MultiQC aggregate report ==========
    log "[QC Step 4/4] Generating MultiQC aggregate report..."
    multiqc "${QC_DIR}" -o "${MULTIQC_DIR}" -n quality_control_report --force || \
        log "WARNING: MultiQC failed or not installed. Skipping aggregate report."
    log "MultiQC report generated (if available)."

    # Create new manifest for cleaned reads
    log "Creating manifest for cleaned sequences..."
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${TEMP_DIR}/cleaned_manifest.tsv"

    for R1 in "${TEMP_DIR}/cleaned/"*_R1.fastq.gz; do
        if [[ -f "${R1}" ]]; then
            SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
            R2="${R1/_R1/_R2}"
            echo -e "${SAMPLE}\t${R1}\t${R2}" >> "${TEMP_DIR}/cleaned_manifest.tsv"
        fi
    done

    # Import cleaned sequences into QIIME2
    log "Importing cleaned sequences into QIIME2..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${TEMP_DIR}/cleaned_manifest.tsv" \
        --output-path "${PROCESSED_DATA_DIR}/fastp-cleaned-seqs.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to import cleaned sequences."

    log "Quality Control completed. Reports in ${QC_DIR}, MultiQC summary in ${MULTIQC_DIR}"
    log "Cleaned sequences: ${PROCESSED_DATA_DIR}/fastp-cleaned-seqs.qza"
}

# This step has been consolidated into step2_quality_control

# Step 3: Remove Host DNA with BBDuk
step3_remove_host() {
    log "Starting Step 3: Removing host DNA with BBDuk"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" "Step 3"; then
        return 0
    fi

    # Check if host reference genome exists
    HOST_REF="${REFERENCE_DB_DIR}/host_genome.fasta"
    if [[ ! -f "${HOST_REF}" ]]; then
        log "WARNING: Host genome reference not found at ${HOST_REF}. Skipping host removal."
        log "To enable host removal, place host reference genome at ${HOST_REF}"
        # Copy input to output if skipping
        cp "${PROCESSED_DATA_DIR}/bbduk-cleaned-seqs.qza" "${PROCESSED_DATA_DIR}/host-removed-seqs.qza"
        return 0
    fi

    # Export sequences from QIIME2 for BBDuk processing
    log "Exporting sequences for host removal..."
    TEMP_DIR="${PROCESSED_DATA_DIR}/temp_host_removal"
    mkdir -p "${TEMP_DIR}"
    
    qiime tools export \
        --input-path "${PROCESSED_DATA_DIR}/bbduk-cleaned-seqs.qza" \
        --output-path "${TEMP_DIR}/exported" || \
        error_exit "Failed to export sequences."

    # Process each sample pair with BBDuk for host removal
    log "Running BBDuk for host DNA removal..."
    mkdir -p "${TEMP_DIR}/host_removed"
    
    for R1 in "${TEMP_DIR}/exported/"*_R1*.fastq*; do
        if [[ -f "${R1}" ]]; then
            R2="${R1/_R1/_R2}"
            SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
            
            log "Removing host DNA from sample: ${SAMPLE}"
            
            bbduk.sh \
                in1="${R1}" \
                in2="${R2}" \
                out1="${TEMP_DIR}/host_removed/${SAMPLE}_R1.fastq.gz" \
                out2="${TEMP_DIR}/host_removed/${SAMPLE}_R2.fastq.gz" \
                ref="${HOST_REF}" \
                k=31 hdist=1 \
                threads=4 \
                stats="${TEMP_DIR}/host_removed/${SAMPLE}_host_stats.txt" || \
                error_exit "BBDuk host removal failed for sample ${SAMPLE}"
        fi
    done

    # Create new manifest
    log "Creating manifest for host-removed sequences..."
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${TEMP_DIR}/host_removed_manifest.tsv"
    
    for R1 in "${TEMP_DIR}/host_removed/"*_R1.fastq.gz; do
        if [[ -f "${R1}" ]]; then
            SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
            R2="${R1/_R1/_R2}"
            echo -e "${SAMPLE}\t${R1}\t${R2}" >> "${TEMP_DIR}/host_removed_manifest.tsv"
        fi
    done

    # Re-import into QIIME2
    log "Re-importing host-removed sequences..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${TEMP_DIR}/host_removed_manifest.tsv" \
        --output-path "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to re-import host-removed sequences."

    # Cleanup
    rm -rf "${TEMP_DIR}"

    log "Host DNA removed successfully with BBDuk. Output: ${PROCESSED_DATA_DIR}/host-removed-seqs.qza"
}

# Step 4: Assemble Contigs
step4_assemble_contigs() {
    log "Starting Step 4: Assembling contigs with ${ASSEMBLER} (${ASSEMBLY_MODE} mode)"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/assembled-contigs.qza" "Step 4"; then
        return 0
    fi

    # Prepare assembly parameters based on mode
    if [[ "${ASSEMBLY_MODE}" == "coassembly" ]]; then
        log "Co-assembly mode: pooling all samples together for assembly"
        COASSEMBLE_FLAG="--p-coassemble"
    else
        log "Per-sample mode: assembling each sample independently"
        COASSEMBLE_FLAG=""
    fi

    # Use host-removed-seqs.qza if exists, otherwise use fastp-cleaned-seqs.qza
    if [[ -f "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" ]]; then
        INPUT_QZA="${PROCESSED_DATA_DIR}/host-removed-seqs.qza"
    else
        log "host-removed-seqs.qza not found, using fastp-cleaned-seqs.qza as input for assembly."
        INPUT_QZA="${PROCESSED_DATA_DIR}/fastp-cleaned-seqs.qza"
    fi

    if [[ "${ASSEMBLER}" == "spades" ]]; then
        # Assembly with metaSPAdes - more accurate but slower
        log "Using metaSPAdes assembler (more accurate, requires more memory)"
        qiime assembly assemble-spades \
            --i-reads "$INPUT_QZA" \
            --p-min-contig-length 500 \
            ${COASSEMBLE_FLAG} \
            --o-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" || \
            error_exit "Failed to assemble contigs with metaSPAdes."
    else
        # Assembly with MEGAHIT - faster and lower memory usage
        log "Using MEGAHIT assembler (faster, lower memory)"
        qiime assembly assemble-megahit \
            --i-reads "$INPUT_QZA" \
            --p-min-contig-len 500 \
            ${COASSEMBLE_FLAG} \
            --o-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" || \
            error_exit "Failed to assemble contigs with MEGAHIT."
    fi

    log "Contigs assembled successfully with ${ASSEMBLER} in ${ASSEMBLY_MODE} mode. Output: ${PROCESSED_DATA_DIR}/assembled-contigs.qza"
    
    # Diagnostic: Export and check contig stats
    log "Checking assembly quality..."
    TEMP_EXPORT="${PROCESSED_DATA_DIR}/temp_contig_check"
    mkdir -p "${TEMP_EXPORT}"
    qiime tools export \
        --input-path "${PROCESSED_DATA_DIR}/assembled-contigs.qza" \
        --output-path "${TEMP_EXPORT}" >/dev/null 2>&1
    
    for contig_file in "${TEMP_EXPORT}"/*.fa "${TEMP_EXPORT}"/*.fna; do
        if [[ -f "${contig_file}" ]]; then
            sample=$(basename "${contig_file}" | sed 's/_contigs.*//')
            total=$(grep -c "^>" "${contig_file}" 2>/dev/null || echo "0")
            if [[ ${total} -gt 0 ]]; then
                avg_len=$(awk '/^>/ {if (seq) print length(seq); seq=""} !/^>/ {seq=seq$0} END {if (seq) print length(seq)}' "${contig_file}" | awk '{sum+=$1; n++} END {if (n>0) print int(sum/n); else print 0}')
                max_len=$(awk '/^>/ {if (seq) print length(seq); seq=""} !/^>/ {seq=seq$0} END {if (seq) print length(seq)}' "${contig_file}" | sort -n | tail -1)
                log "${sample}: ${total} contigs, avg=${avg_len}bp, max=${max_len}bp"
            fi
        fi
    done
    rm -rf "${TEMP_EXPORT}"
}

# Step 4b: Map reads to contigs (required for binning)
step4b_map_reads() {
    log "Starting Step 4b: Mapping reads to contigs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mapped-reads.qza" "Step 4b"; then
        return 0
    fi

    # Index contigs
    qiime assembly index-contigs \
        --i-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" \
        --o-index "${PROCESSED_DATA_DIR}/contigs-index.qza" || \
        error_exit "Failed to index contigs."

    # Map reads to contigs
    qiime assembly map-reads \
        --i-index "${PROCESSED_DATA_DIR}/contigs-index.qza" \
        --i-reads "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" \
        --o-alignment-maps "${PROCESSED_DATA_DIR}/mapped-reads.qza" || \
        error_exit "Failed to map reads to contigs."

    log "Reads mapped successfully. Output: ${PROCESSED_DATA_DIR}/mapped-reads.qza"
}

# Step 5: Bin Contigs (per-sample binning with MetaBAT2)
step5_bin_contigs() {
    log "Starting Step 5: Per-sample binning with MetaBAT2"

    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/binned-contigs.qza" "Step 5"; then
        return 0
    fi

    # Map reads to contigs (required for MetaBAT2)
    log "Mapping reads to contigs for MetaBAT2..."
    qiime assembly map-reads \
        --i-index "${PROCESSED_DATA_DIR}/contigs-index.qza" \
        --i-reads "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" \
        --o-alignment-maps "${PROCESSED_DATA_DIR}/mapped-reads.qza" || \
        error_exit "Failed to map reads to contigs."

    # Run MetaBAT2 for binning
    log "Running MetaBAT2 for binning..."
    qiime metagenomics bin-contigs-metabat \
        --i-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" \
        --i-alignment-maps "${PROCESSED_DATA_DIR}/mapped-reads.qza" \
        --o-bins "${PROCESSED_DATA_DIR}/binned-contigs.qza" || \
        error_exit "Failed to bin contigs with MetaBAT2."

    log "Per-sample binning completed successfully with MetaBAT2. Output: ${PROCESSED_DATA_DIR}/binned-contigs.qza"
}

# Step 5b: Evaluate MAGs
step5b_evaluate_mags() {
    log "Starting Step 5b: Evaluating MAGs with CheckM2"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mags-evaluation.qza" "Step 5b"; then
        return 0
    fi

    qiime quality-control evaluate-checkm2 \
        --i-bins "${PROCESSED_DATA_DIR}/binned-contigs.qza" \
        --o-results "${PROCESSED_DATA_DIR}/mags-evaluation.qza" \
        --o-visualization "${PROCESSED_DATA_DIR}/mags-evaluation.qzv" || \
        error_exit "Failed to evaluate MAGs with CheckM2."

    log "MAGs evaluated successfully with CheckM2. Outputs: ${PROCESSED_DATA_DIR}/mags-evaluation.qza, ${PROCESSED_DATA_DIR}/mags-evaluation.qzv"
}

# Step 6: Classify MAG Taxonomy with GTDB-Tk
step6_classify_taxonomy() {
    log "Starting Step 6: Taxonomic classification of MAGs with GTDB-Tk"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mags-taxonomy.qza" "Step 6"; then
        return 0
    fi

    # Check if GTDB database exists
    GTDB_DATA="${REFERENCE_DB_DIR}/gtdb"
    if [[ ! -d "${GTDB_DATA}" ]]; then
        log "WARNING: GTDB database not found at ${GTDB_DATA}."
        log "To enable GTDB-Tk classification, download GTDB database:"
        log "  wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz"
        log "  tar -xzf gtdbtk_data.tar.gz -C ${REFERENCE_DB_DIR}"
        log "Skipping taxonomic classification."
        return 0
    fi

    log "Running GTDB-Tk for MAG classification..."
    qiime moshpit classify-gtdbtk \
        --i-bins "${PROCESSED_DATA_DIR}/binned-contigs.qza" \
        --p-gtdbtk-data-path "${GTDB_DATA}" \
        --p-num-cpus 4 \
        --o-classification "${PROCESSED_DATA_DIR}/mags-taxonomy.qza" \
        --o-summary "${PROCESSED_DATA_DIR}/mags-taxonomy-summary.qzv" || \
        error_exit "Failed to classify MAGs with GTDB-Tk."

    log "MAG taxonomy classification completed. Outputs: ${PROCESSED_DATA_DIR}/mags-taxonomy.qza, ${PROCESSED_DATA_DIR}/mags-taxonomy-summary.qzv"
}

# Step 7a: Predict genes from MAGs
step7a_predict_genes() {
    log "Starting Step 7a: Predicting genes"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/predicted-genes.qza" "Step 7a"; then
        return 0
    fi

    qiime annotate predict-genes-prodigal \
        --i-input "${PROCESSED_DATA_DIR}/binned-contigs.qza" \
        --o-gene-sequences "${PROCESSED_DATA_DIR}/predicted-genes.qza" || \
        error_exit "Failed to predict genes."

    log "Genes predicted successfully. Output: ${PROCESSED_DATA_DIR}/predicted-genes.qza"
}


# Step 7c: Annotate MAGs with eggNOG
step7c_annotate_mags() {
    log "Starting Step 7c: Annotating MAGs with eggNOG"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mags-annotations.qza" "Step 7c"; then
        return 0
    fi

    qiime annotate map-eggnog \
        --i-ortholog-annotations "${PROCESSED_DATA_DIR}/orthologs.qza" \
        --o-eggnog-annotations "${PROCESSED_DATA_DIR}/mags-annotations.qza" || \
        error_exit "Fail vbed to annotate MAGs."

    log "MAGs annotated successfully. Output: ${PROCESSED_DATA_DIR}/mags-annotations.qza"
}

# Parse arguments
TEST_MODE=false
START_STEP=0
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test)
            TEST_MODE=true
            shift
            ;;
        --start-step|--start)
            if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
                START_STEP="$2"
                shift 2
            else
                log "ERROR: --start-step requires a numeric argument"
                exit 1
            fi
            ;;
        --assembler)
            if [[ -n "$2" && "$2" =~ ^(megahit|spades)$ ]]; then
                ASSEMBLER="$2"
                shift 2
            else
                log "ERROR: --assembler requires 'megahit' or 'spades' as argument"
                exit 1
            fi
            ;;
        --assembly-mode)
            if [[ -n "$2" && "$2" =~ ^(per-sample|coassembly)$ ]]; then
                ASSEMBLY_MODE="$2"
                shift 2
            else
                log "ERROR: --assembly-mode requires 'per-sample' or 'coassembly' as argument"
                exit 1
            fi
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --test                  Run in test mode"
            echo "  --start-step N          Start from step N"
            echo "  --assembler TYPE        Choose assembler: megahit (default, fast) or spades (accurate)"
            echo "  --assembly-mode MODE    Choose mode: per-sample (default) or coassembly (pool all samples)"
            echo "  -h, --help              Show this help message"
            echo ""
            echo "Assembly modes:"
            echo "  per-sample   - Assemble each sample independently (preserves sample info)"
            echo "  coassembly   - Pool all samples together (better for low-abundance organisms)"
            exit 0
            ;;
        *)
            log "Unknown option: $1"
            log "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Adjust directories for test mode
if $TEST_MODE; then
    log "Running in test mode. Using test dataset."
    RAW_DATA_DIR="${DATA_DIR}/raw_data/shotgun_test"
    PROCESSED_DATA_DIR="${DATA_DIR}/processed_data_test"
    mkdir -p "${RAW_DATA_DIR}" "${PROCESSED_DATA_DIR}"
fi

# Helper to run a step only if its step number is >= START_STEP
run_step() {
    local step_num="$1"
    local func_name="$2"
    if (( step_num < START_STEP )); then
        log "Skipping ${func_name} (step ${step_num}) because START_STEP=${START_STEP}"
        return 0
    fi
    ${func_name}
}

# Main Pipeline Execution (Updated)
main() {
    log "Starting QIIME2 Shotgun Data Analysis Pipeline (START_STEP=${START_STEP})"

    run_step 1 step1_import_data
    run_step 2 step2_quality_control
    run_step 3 step3_remove_host
    run_step 4 step4_assemble_contigs
    run_step 5 step4b_map_reads
    run_step 6 step5_bin_contigs
    run_step 7 step5b_evaluate_mags
    run_step 8 step6_classify_taxonomy
    run_step 9 step7a_predict_genes
    # run_step 10 step7b_search_orthologs (removed)
    run_step 11 step7c_annotate_mags

    log "Pipeline completed successfully."
}

main "$@"
