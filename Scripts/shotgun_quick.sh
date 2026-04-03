#!/bin/bash
set -euo pipefail

# ============================================================================
# Quick Metagenomic Search - Methane Gene Detection Pipeline
# ============================================================================

# --- Defaults ---
THREADS=${THREADS:-4}
MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-1}
RAW_DATA_SUBDIR=${RAW_DATA_SUBDIR:-}

# --- Resolve project paths from this script's location ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"

# --- Temp dir ---
TMPDIR_BASE="${TMPDIR_BASE:-${BASE_DIR}/Data/temp}"
mkdir -p "$TMPDIR_BASE"
export TMPDIR="$TMPDIR_BASE" TEMP="$TMPDIR_BASE" TMP="$TMPDIR_BASE"

# ============================================================================
# Conda environment check
# ============================================================================
if ! command -v conda &>/dev/null; then
    source "${HOME}/.bashrc" 2>/dev/null || true
    CONDA_BASE="${CONDA_BASE:-$(conda info --base 2>/dev/null)}"
    source "${CONDA_BASE}/etc/profile.d/conda.sh" 2>/dev/null || true
fi

if [[ "$CONDA_DEFAULT_ENV" != "quick_search" ]]; then
    if [[ -n "${SLURM_JOB_ID:-}${PBS_JOBID:-}${PBS_JOBNAME:-}" ]]; then
        echo "HPC job detected. Activating quick_search..."
        conda activate quick_search || { echo "ERROR: Failed to activate quick_search"; exit 1; }
    else
        echo "ERROR: Activate the 'quick_search' conda environment first:"
        echo "  conda activate quick_search"
        exit 1
    fi
fi

# ============================================================================
# Argument parsing
# ============================================================================
SAMPLE_LIST=()
SAMPLE_FILE=""
AUTO_DETECT=false
FORCE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        step0|--setup)
            echo "1. conda create -n quick_search -c bioconda -c conda-forge python=3.10 kraken2 fastp diamond bowtie2"
            echo "2. conda activate quick_search"
            echo "3. kraken2-build --standard --threads 4 --db Data/reference_dbs/kraken2_db"
            echo "4. Place DIAMOND DB at: Data/reference_dbs/DIAMOND/methane_master_db.dmnd"
            exit 0 ;;
        --sample-list) SAMPLE_FILE="$2"; shift 2 ;;
        --auto)        AUTO_DETECT=true; shift ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --parallel)    MAX_PARALLEL_SAMPLES="$2"; shift 2 ;;
        --base-dir)    BASE_DIR="$2"; shift 2 ;;
        --force)       FORCE=true; shift ;;
        --help|-h)
            cat <<'EOF'
USAGE: bash Scripts/shotgun_quick.sh [OPTIONS] [SAMPLE_IDs...]

OPTIONS:
  step0, --setup        Show setup instructions
  --sample-list FILE    Read sample IDs from file (one per line)
  --auto                Auto-detect samples in raw_data directory
  --threads N           Threads (default: 4)
  --parallel N          Parallel samples (default: 1)
  --base-dir PATH       Base directory (default: auto-detect)
  --force               Re-run all steps even if outputs exist
  --help, -h            This help

EXAMPLES:
  bash Scripts/shotgun_quick.sh 53394
  bash Scripts/shotgun_quick.sh --auto --threads 16
  bash Scripts/shotgun_quick.sh --sample-list Config/samples.txt

REQUIREMENTS:
  - Conda env: quick_search (fastp, diamond, orfm, hmmer, kraken2)
  - Kraken2 DB:  Data/reference_dbs/kraken2_db/
  - DIAMOND DB:  Data/reference_dbs/DIAMOND/methane_master_db.dmnd
  - Input:       Data/raw_data/<SAMPLE>_R1*.fastq.gz
EOF
            exit 0 ;;
        -*) echo "ERROR: Unknown option: $1. Use --help."; exit 1 ;;
        *)  SAMPLE_LIST+=("$1"); shift ;;
    esac
done

# ============================================================================
# Base dir is resolved from the script location unless --base-dir is provided
# ============================================================================

# ============================================================================
# Collect samples
# ============================================================================
RAW_DATA_DIR="${BASE_DIR}/Data/raw_data/${RAW_DATA_SUBDIR}"

if $AUTO_DETECT; then
    echo "Auto-detecting samples in ${RAW_DATA_DIR}..."
    while IFS= read -r f; do
        SAMPLE_LIST+=("$(basename "$f" | sed -E 's/_R[12].*.fastq.gz$//')")
    done < <(find "$RAW_DATA_DIR" -name "*_R1*.fastq.gz" -type f)
elif [ -n "$SAMPLE_FILE" ]; then
    [ -f "$SAMPLE_FILE" ] || { echo "ERROR: Sample file not found: $SAMPLE_FILE"; exit 1; }
    while IFS= read -r line || [ -n "$line" ]; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        SAMPLE_LIST+=("$line")
    done < "$SAMPLE_FILE"
fi

[ ${#SAMPLE_LIST[@]} -gt 0 ] || {
    echo "ERROR: No samples provided. Use --help for usage."; exit 1
}

# ============================================================================
# Prerequisites check
# ============================================================================
echo "Checking tools..."
MISSING=()
for tool in fastp diamond orfm hmmsearch bc; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
if [ ${#MISSING[@]} -gt 0 ]; then
    echo "ERROR: Missing tools: ${MISSING[*]}"
    echo "  conda install -c bioconda -c conda-forge fastp diamond orfm hmmer"
    exit 1
fi
echo "  ✓ fastp  $(fastp --version 2>&1 | head -1)"
echo "  ✓ diamond  $(diamond --version 2>&1 | head -1)"

# ============================================================================
# Paths and database check
# ============================================================================
FASTP_DIR="${BASE_DIR}/Data/processed_data/fastp_cleaned"
HOST_REMOVED_DIR="${BASE_DIR}/Data/processed_data/host_removed"
KRAKEN2_DIR="${BASE_DIR}/Results/kraken2_output"
DIAMOND_DIR="${BASE_DIR}/Results/functional_analysis"
LOG_DIR="${BASE_DIR}/Logs"
mkdir -p "$FASTP_DIR" "$HOST_REMOVED_DIR" "$KRAKEN2_DIR" "$DIAMOND_DIR" "$LOG_DIR"

# Auto-detect Kraken2 DB
if [ -z "${KRAKEN2_DB:-}" ]; then
    if   [ -f "${BASE_DIR}/Data/reference_dbs/kraken2_db/hash.k2d" ]; then
        KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs/kraken2_db"
    elif [ -f "${BASE_DIR}/Data/reference_dbs/hash.k2d" ]; then
        KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs"
    else
        KRAKEN2_DB="${BASE_DIR}/Data/reference_dbs/kraken2_db"
    fi
fi
METHANE_DB="${METHANE_DB:-${BASE_DIR}/Data/reference_dbs/DIAMOND/methane_master_db}"
BT2_INDEX="${BASE_DIR}/Data/reference_dbs/host_genomes/plant_host"

echo "Checking databases..."
DB_MISSING=false
[ -f "${KRAKEN2_DB}/hash.k2d" ] && echo "  ✓ Kraken2 DB: ${KRAKEN2_DB}" \
    || { echo "  Kraken2 DB NOT found at: ${KRAKEN2_DB}"; DB_MISSING=true; }
{ [ -f "${METHANE_DB}.dmnd" ] || [ -f "${METHANE_DB}" ]; } && echo "  ✓ DIAMOND DB: ${METHANE_DB}" \
    || { echo "  DIAMOND DB NOT found at: ${METHANE_DB}.dmnd"; DB_MISSING=true; }
$DB_MISSING && { echo "ERROR: Required databases missing. See --help."; exit 1; }

echo "Ready to process ${#SAMPLE_LIST[@]} sample(s) | threads=${THREADS} | base=${BASE_DIR} | force=${FORCE}"
echo ""

# ============================================================================
# Checkpoint helper
# skip_if_done FILE [FILE2...] → returns 0 (skip) if all files exist and FORCE=false
# ============================================================================
skip_if_done() {
    $FORCE && return 1   # --force: never skip
    for f in "$@"; do
        [ -f "$f" ] || return 1   # any file missing → don't skip
    done
    return 0   # all files exist → skip
}

# ============================================================================
# Process single sample
# ============================================================================
process_sample() {
    local SAMPLE_ID=$1
    local START_TIME LOG_FILE
    START_TIME=$(date +%s)
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_quick_search.log"

    (
    echo "========================================================"
    echo "Processing: ${SAMPLE_ID}  [$(date)]"
    echo "========================================================"

    # --- Locate input files ---
    if [ -f "${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz" ]; then
        R1_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz"
        R2_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    else
        R1_RAW=$(find "$RAW_DATA_DIR" \( -name "${SAMPLE_ID}_R1_*.fastq.gz" -o -name "${SAMPLE_ID}_R1.fastq.gz" \) | head -1)
        R2_RAW=$(find "$RAW_DATA_DIR" \( -name "${SAMPLE_ID}_R2_*.fastq.gz" -o -name "${SAMPLE_ID}_R2.fastq.gz" \) | head -1)
    fi
    [ -f "${R1_RAW:-}" ] && [ -f "${R2_RAW:-}" ] \
        || { echo "ERROR: Input files not found for ${SAMPLE_ID}"; return 1; }
    echo "Input: $(basename "$R1_RAW") / $(basename "$R2_RAW")"

    # --- Output paths ---
    R1_CLEAN="${FASTP_DIR}/${SAMPLE_ID}_R1_cleaned.fastq.gz"
    R2_CLEAN="${FASTP_DIR}/${SAMPLE_ID}_R2_cleaned.fastq.gz"
    FASTP_JSON="${FASTP_DIR}/${SAMPLE_ID}_fastp.json"
    FASTP_HTML="${FASTP_DIR}/${SAMPLE_ID}_fastp.html"
    KRAKEN2_OUTPUT="${KRAKEN2_DIR}/kraken2_output_${SAMPLE_ID}"
    DIAMOND_OUTPUT="${DIAMOND_DIR}/${SAMPLE_ID}_combined_methane_hits.txt"
    DIAMOND_SUMMARY="${DIAMOND_DIR}/${SAMPLE_ID}_methane_summary.txt"
    FILTERED_HITS="${DIAMOND_DIR}/${SAMPLE_ID}_filtered_methane_hits.txt"
    GENE_COUNTS="${DIAMOND_DIR}/${SAMPLE_ID}_gene_counts.txt"
    GENE_RPKM="${DIAMOND_DIR}/${SAMPLE_ID}_gene_rpkm.txt"

    # --------------------------------------------------------
    # STEP 1: fastp Quality Control
    # --------------------------------------------------------
    echo "--- STEP 1: fastp QC ---"
    if skip_if_done "$R1_CLEAN" "$R2_CLEAN"; then
        echo "SKIP: fastp outputs already exist"
    else
        fastp -i "$R1_RAW" -I "$R2_RAW" -o "$R1_CLEAN" -O "$R2_CLEAN" \
            --dedup --dup_calc_accuracy 3 --detect_adapter_for_pe \
            --qualified_quality_phred 20 --length_required 50 \
            --thread "$THREADS" --json "$FASTP_JSON" --html "$FASTP_HTML" \
            2>&1 | tee "${LOG_DIR}/${SAMPLE_ID}_fastp.log"
        echo "fastp done. QC report: ${FASTP_HTML}"
    fi

    # --------------------------------------------------------
    # --------------------------------------------------------
    # STEP 2: Host Removal (Bowtie2)
    # --------------------------------------------------------
    echo "--- STEP 2: Host Removal ---"
    R1_FOR_ANALYSIS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R1_host_removed.fastq.gz"
    R2_FOR_ANALYSIS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_R2_host_removed.fastq.gz"
    HOST_STATS="${HOST_REMOVED_DIR}/${SAMPLE_ID}_host_removal_stats.txt"
    
    if [ -f "${BT2_INDEX}.1.bt2" ]; then
        if skip_if_done "$R1_FOR_ANALYSIS" "$R2_FOR_ANALYSIS"; then
            echo "SKIP: Host-removed files already exist"
        else
            echo "Removing host reads using Bowtie2 index: ${BT2_INDEX}"
            bowtie2 -x "$BT2_INDEX" -1 "$R1_CLEAN" -2 "$R2_CLEAN" \
                --un-conc-gz "${HOST_REMOVED_DIR}/${SAMPLE_ID}_R%_host_removed.fastq.gz" \
                --threads "$THREADS" --very-sensitive-local -S /dev/null \
                2>&1 | tee "$HOST_STATS"
            echo "Host removal complete. Stats: ${HOST_STATS}"
        fi
    else
        echo "WARNING: Bowtie2 index not found at ${BT2_INDEX}"
        echo "Symlinking cleaned files as host-removed (no actual removal)"
        ln -sf "$R1_CLEAN" "$R1_FOR_ANALYSIS"
        ln -sf "$R2_CLEAN" "$R2_FOR_ANALYSIS"
        echo "No host removal performed" > "$HOST_STATS"
    fi
    # --------------------------------------------------------
    # STEP 3: Kraken2 Taxonomic Profiling
    # --------------------------------------------------------
    echo "--- STEP 3: Kraken2 ---"
    # Find kraken2 binary
    if command -v kraken2 &>/dev/null; then
        KRAKEN2_BIN="kraken2"
    else
        KRAKEN2_BIN="${CONDA_BASE}/envs/kraken_env/bin/kraken2"
        [ -f "$KRAKEN2_BIN" ] || { echo "ERROR: kraken2 not found"; return 1; }
        echo "WARNING: Using kraken2 from kraken_env (consider installing in quick_search)"
    fi
    $KRAKEN2_BIN --version &>/dev/null || { echo "ERROR: kraken2 not executable"; return 1; }

    SAMPLE_TMPDIR="${TMPDIR}/kraken2_${SAMPLE_ID}_$$"
    mkdir -p "$SAMPLE_TMPDIR"
    export TMPDIR="$SAMPLE_TMPDIR"
    export KRAKEN2_DB_PATH="$KRAKEN2_DB"
    KRAKEN2_THREADS=${KRAKEN2_THREADS:-$THREADS}
    export OMP_NUM_THREADS=1
    ulimit -s 65536 2>/dev/null || true

    mkdir -p "$KRAKEN2_OUTPUT"
    if skip_if_done "${KRAKEN2_OUTPUT}/kraken2_report.txt"; then
        echo "SKIP: Kraken2 output already exists"
    else
        $KRAKEN2_BIN --db "$KRAKEN2_DB" --paired "$R1_FOR_ANALYSIS" "$R2_FOR_ANALYSIS" \
            --output "${KRAKEN2_OUTPUT}/kraken2_output.txt" \
            --report "${KRAKEN2_OUTPUT}/kraken2_report.txt" \
            --threads "$KRAKEN2_THREADS" --use-names \
            || { echo "ERROR: Kraken2 failed"; return 1; }
        echo "Kraken2 done. Report: ${KRAKEN2_OUTPUT}/kraken2_report.txt"
    fi
    rm -rf "$SAMPLE_TMPDIR" 2>/dev/null || true
    export TMPDIR="$TMPDIR_BASE"

    # --------------------------------------------------------
    # --------------------------------------------------------
    # STEP 4: Functional Gene Analysis (DIAMOND + Stats + RPKM)
    # --------------------------------------------------------
    echo "--- STEP 4: Functional Gene Analysis ---"
    
    # Check if all final outputs exist (skip entire analysis if complete)
    if skip_if_done "$DIAMOND_OUTPUT" "$DIAMOND_SUMMARY" "$GENE_RPKM"; then
        echo "SKIP: All functional analysis outputs already exist"
        # Load existing stats for final report
        TOTAL_HITS=$(wc -l < "$DIAMOND_OUTPUT")
        MMO_HITS=$(grep -ic "methane monooxygenase" "$DIAMOND_OUTPUT" || true)
        PMO_HITS=$(grep -iEc "PmoA|pMMO" "$DIAMOND_OUTPUT" || true)
        TOTAL_METHANE=$((MMO_HITS + PMO_HITS))
        FILTERED_COUNT=$(wc -l < "$FILTERED_HITS")
        TOTAL_READS=$(( $(zcat "$R1_FOR_ANALYSIS" "$R2_FOR_ANALYSIS" | wc -l) / 4 ))
        TOTAL_READS_M=$(awk "BEGIN {printf \"%.4f\", ${TOTAL_READS}/1000000}")
        TOTAL_MAPPED_READS=$(awk '{print $1}' "$DIAMOND_OUTPUT" | sort -u | wc -l)
    else
        # 4.1: DIAMOND alignment
        echo "  4.1: Running DIAMOND alignment..."
        DIAMOND_TMPDIR="${TMPDIR}/diamond_${SAMPLE_ID}_$$"
        mkdir -p "$DIAMOND_TMPDIR"
        gunzip -c "$R1_FOR_ANALYSIS" "$R2_FOR_ANALYSIS" | \
            DIAMOND_TMPDIR="$DIAMOND_TMPDIR" diamond blastx \
            -d "$METHANE_DB" -q - -o "$DIAMOND_OUTPUT" \
            --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
            --threads "$THREADS" --sensitive \
            --query-cover 80 --min-score 40 --max-target-seqs 5 --evalue 1e-5 \
            || { echo "ERROR: DIAMOND failed"; return 1; }
        rm -rf "$DIAMOND_TMPDIR" 2>/dev/null || true
        echo "  DIAMOND alignment complete: ${DIAMOND_OUTPUT}"
        
        # 4.2: Generate summary statistics
        echo "  4.2: Generating summary statistics..."
        TOTAL_HITS=$(wc -l < "$DIAMOND_OUTPUT")
        MMO_HITS=$(grep -ic "methane monooxygenase" "$DIAMOND_OUTPUT" || true)
        PMO_HITS=$(grep -iEc "PmoA|pMMO" "$DIAMOND_OUTPUT" || true)
        TOTAL_METHANE=$((MMO_HITS + PMO_HITS))
        
        cat > "$DIAMOND_SUMMARY" <<REPORT
METHANE METABOLISM GENE SEARCH RESULTS - ${SAMPLE_ID} [$(date +%Y-%m-%d)]
DB: methane_master_db | query-cover >80% | bitscore >40 | E-value <1e-5

Total alignments:            ${TOTAL_HITS}
Methane monooxygenase (MMO): ${MMO_HITS}
PmoA family (pMMO):          ${PMO_HITS}
Total methanotroph hits:     ${TOTAL_METHANE}
REPORT
        
        if [ "$TOTAL_METHANE" -gt 0 ]; then
            cat >> "$DIAMOND_SUMMARY" <<INTERP
METHANOTROPHS DETECTED: Aerobic methane-oxidizing bacteria present.
  pMMO = primary methane oxidation; MMO = alternative pathway.
INTERP
        else
            echo "NO SIGNIFICANT METHANE METABOLISM GENES DETECTED." >> "$DIAMOND_SUMMARY"
        fi
        echo "  Summary report: ${DIAMOND_SUMMARY}"
        
        # 4.3: Calculate gene abundance (RPKM)
        echo "  4.3: Calculating gene abundance (RPKM)..."
        grep -iE "methane|pmo|mmo" "$DIAMOND_OUTPUT" \
            | grep -viE "ribosom|ribokinase|ribonucleoside|ammonia" > "$FILTERED_HITS"
        FILTERED_COUNT=$(wc -l < "$FILTERED_HITS")
        echo "  Filtered hits: ${FILTERED_COUNT}"
        
        awk '{print $7}' "$FILTERED_HITS" | sort | uniq -c | sort -rn > "$GENE_COUNTS"
        echo "  Top 10 genes:"; head -10 "$GENE_COUNTS"
        
        TOTAL_MAPPED_READS=$(awk '{print $1}' "$DIAMOND_OUTPUT" | sort -u | wc -l)
        TOTAL_READS=$(( $(zcat "$R1_FOR_ANALYSIS" "$R2_FOR_ANALYSIS" | wc -l) / 4 ))
        TOTAL_READS_M=$(awk "BEGIN {printf \"%.4f\", ${TOTAL_READS}/1000000}")
        echo "  Total reads: ${TOTAL_READS} (${TOTAL_READS_M}M) | Mapped: ${TOTAL_MAPPED_READS}"
        
        awk -v total_reads="$TOTAL_READS_M" -v filt="$FILTERED_HITS" '
        BEGIN { print "Gene_Name\tRead_Count\tAvg_Length\tRPKM" }
        {
            gsub(/^[[:space:]]+/, ""); split($0, a, " ")
            count = a[1]; gene = ""
            for (i=2; i<=length(a); i++) gene = gene " " a[i]
            gsub(/^[[:space:]]+/, "", gene)
            cmd = "grep -F \"" gene "\" " filt " | awk '"'"'{s+=$4;n++} END{print (n>0)?s/n:1000}'"'"'"
            cmd | getline avg_len; close(cmd)
            printf "%s\t%d\t%.1f\t%.4f\n", gene, count, avg_len, (count*1000)/(avg_len*total_reads)
        }' "$GENE_COUNTS" > "$GENE_RPKM"
        echo "  Top 10 by RPKM:"; head -11 "$GENE_RPKM" | tail -10
        echo "  RPKM analysis complete: ${GENE_RPKM}"
    fi
    # STEP 5: Results
    # --------------------------------------------------------
    local END_TIME ELAPSED
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo ""
    echo "=== COMPLETE: ${SAMPLE_ID} | $(( ELAPSED/60 ))m $(( ELAPSED%60 ))s ==="
    echo "  QC:       ${FASTP_HTML}"
    echo "  Kraken2:  ${KRAKEN2_OUTPUT}/"
    echo "  DIAMOND:  ${DIAMOND_OUTPUT}"
    echo "  RPKM:     ${GENE_RPKM}"
    echo "  Hits: total=${TOTAL_HITS} | filtered=${FILTERED_COUNT} | MMO=${MMO_HITS} | pMMO=${PMO_HITS}"

    ) > >(tee -a "$LOG_FILE") 2>&1
    return $?
}

# ============================================================================
# Main execution loop
# ============================================================================
MAIN_START=$(date +%s)
SUCCESS_COUNT=0 FAIL_COUNT=0 FAILED_SAMPLES=()

if [ "$MAX_PARALLEL_SAMPLES" -gt 1 ]; then
    if command -v parallel &>/dev/null; then
        echo "Parallel mode: ${MAX_PARALLEL_SAMPLES} jobs..."
        export -f process_sample skip_if_done
        export BASE_DIR THREADS RAW_DATA_DIR FASTP_DIR KRAKEN2_DIR DIAMOND_DIR LOG_DIR \
               KRAKEN2_DB METHANE_DB KRAKEN2_THREADS RAW_DATA_SUBDIR TMPDIR_BASE BT2_INDEX FORCE
        printf '%s\n' "${SAMPLE_LIST[@]}" | \
            parallel -j "$MAX_PARALLEL_SAMPLES" --line-buffer \
            'if process_sample {}; then echo "✓ {} done"; else echo "✗ {} failed"; fi'
        for S in "${SAMPLE_LIST[@]}"; do
            grep -q "COMPLETE:" "${LOG_DIR}/${S}_quick_search.log" 2>/dev/null \
                && { (( SUCCESS_COUNT++ )) || true; } || { (( FAIL_COUNT++ )) || true; FAILED_SAMPLES+=("$S"); }
        done
    else
        echo "WARNING: GNU parallel not found. Falling back to sequential."
        MAX_PARALLEL_SAMPLES=1
    fi
fi

if [ "$MAX_PARALLEL_SAMPLES" -eq 1 ]; then
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo ""
        echo "### Processing ${SAMPLE_ID} ($(( SUCCESS_COUNT+FAIL_COUNT+1 ))/${#SAMPLE_LIST[@]}) ###"
        if process_sample "$SAMPLE_ID"; then
            (( SUCCESS_COUNT++ )) || true
            echo "✓ ${SAMPLE_ID} done"
        else
            (( FAIL_COUNT++ )) || true
            FAILED_SAMPLES+=("$SAMPLE_ID")
            echo "✗ ${SAMPLE_ID} failed"
        fi
    done
fi

# ============================================================================
# Final summary
# ============================================================================
TOTAL_ELAPSED=$(( $(date +%s) - MAIN_START ))
echo ""
echo "========================================================"
echo "PIPELINE SUMMARY"
echo "  Samples:    ${#SAMPLE_LIST[@]}"
echo "  Successful: ${SUCCESS_COUNT}"
echo "  Failed:     ${FAIL_COUNT}"
echo "  Time:       $(( TOTAL_ELAPSED/3600 ))h $(( TOTAL_ELAPSED%3600/60 ))m $(( TOTAL_ELAPSED%60 ))s"
echo "========================================================"

if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "Failed samples: ${FAILED_SAMPLES[*]}"
    echo "Check logs in: ${LOG_DIR}/"
    exit 1
fi

echo "All done! Outputs:"
echo "  fastp:    Data/processed_data/fastp_cleaned/"
echo "  Kraken2:  Data/processed_data/kraken2_output/"
echo "  Methane:  Data/functional_analysis/methane_genes/"
echo "  Logs:     Logs/"
