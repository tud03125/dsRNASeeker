#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PYTHON_EXE="${PYTHON_EXE:-python3}"
RSCRIPT_EXE="${RSCRIPT_EXE:-Rscript}"

MATRIX="${MATRIX:?Set MATRIX to the GSE162876 gene-count matrix}"

GTF="${GTF:?Set GTF to the mm39 annotation GTF}"

DSRROOT="${DSRROOT:?Set DSRROOT to the GSE162876 dsRNASeeker output root}"
OUTROOT="${OUTROOT:-$APP_DIR/reproduced}"

PREP="$SCRIPT_DIR/prepare_salmon_gene_counts.py"
GSEA="$SCRIPT_DIR/run_gene_gsea_candidate_link.R"

mkdir -p "$OUTROOT"

# ============================================================
# PREFLIGHT
# ============================================================

for F in \
    "$MATRIX" \
    "$GTF" \
    "$PREP" \
    "$GSEA"
do
    [[ -s "$F" ]] || {
        echo "[ERROR] Missing required file: $F"
        exit 1
    }
done

echo "===== VERIFY MOUSE HALLMARK FIX ====="

grep -nE \
    'hallmark_collection|db_species.*MM|collection.*MH' \
    "$GSEA" || {
        echo "[ERROR] Could not verify mouse Hallmark collection handling"
        exit 1
    }

# ============================================================
# FUNCTION
# ============================================================

run_one() {

    local CASE="$1"
    local DATASET="$2"
    local SUMMARY="$3"
    local CASE_REGEX="$4"

    local OUTDIR="$OUTROOT/$DATASET"

    mkdir -p "$OUTDIR"

    [[ -s "$SUMMARY" ]] || {
        echo "[ERROR] Missing dsRNASeeker summary: $SUMMARY"
        exit 1
    }

    echo
    echo "============================================================"
    echo "[INFO] $DATASET"
    echo "[INFO] Case: $CASE"
    echo "[INFO] Control: FPC_Control"
    echo "[INFO] Summary: $SUMMARY"
    echo "============================================================"

    "$PYTHON_EXE" "$PREP" \
        --matrix "$MATRIX" \
        --gtf "$GTF" \
        --outdir "$OUTDIR" \
        --case-label "$CASE" \
        --control-label FPC_Control \
        --case-regex "$CASE_REGEX" \
        --control-regex '^RNA-seq_IP_FPC_Control_'

    "$RSCRIPT_EXE" "$GSEA" \
        --counts "$OUTDIR/gene_raw_counts.tsv" \
        --samples "$OUTDIR/gene_samples.tsv" \
        --summary "$SUMMARY" \
        --case "$CASE" \
        --control FPC_Control \
        --species mouse \
        --dataset "$DATASET" \
        --outdir "$OUTDIR" \
        --min-count 10 \
        --min-size 15 \
        --max-size 500 \
        --top-n 15
}

# ============================================================
# FPC MILD
# ============================================================

run_one \
    FPC_Mild \
    GSE162876_FPC_Mild_vs_Control \
    "$DSRROOT/FPC_Mild_vs_Control/inverted/summary/TEpair_dsRNA_master.summary.with_RI.csv" \
    '^RNA-seq_IP_FPC_Mild_'

# ============================================================
# FPC ADVANCED
# ============================================================

run_one \
    FPC_Advanced \
    GSE162876_FPC_Advanced_vs_Control \
    "$DSRROOT/FPC_Advanced_vs_Control/inverted/summary/TEpair_dsRNA_master.summary.with_RI.csv" \
    '^RNA-seq_IP_FPC_Advanced_'

echo
echo "============================================================"
echo "[PASS] Both candidate-linked FPC GSEA analyses completed"
echo "[INFO] Outputs: $OUTROOT"
echo "============================================================"
