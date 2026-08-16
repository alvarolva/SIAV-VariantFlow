#!/usr/bin/env bash
set -euo pipefail

# Quality control and trimming for paired-end Illumina reads
# H1N2 influenza study
#
# Input:
#   data/fastq/<sample>_R1.fastq.gz
#   data/fastq/<sample>_R2.fastq.gz
#
# Output:
#   data/trimmed/<sample>_R1.fastq.gz
#   data/trimmed/<sample>_R2.fastq.gz
#   data/trimmed/unpaired/<sample>_R1.unpaired.fastq.gz
#   data/trimmed/unpaired/<sample>_R2.unpaired.fastq.gz
#
# QC output:
#   results/qc/raw_fastqc/
#   results/qc/raw_multiqc/
#   results/qc/trimmed_fastqc/
#   results/qc/trimmed_multiqc/
#
# Requirements:
#   trimmomatic
#   fastqc
#   multiqc
#
# Usage:
#   bash scripts/01_quality_control.sh
#
# Optional:
#   THREADS=8 bash scripts/01_quality_control.sh


THREADS="${THREADS:-4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RAW_DIR="${REPO_ROOT}/data/fastq"

TRIMMED_DIR="${REPO_ROOT}/data/trimmed"
UNPAIRED_DIR="${TRIMMED_DIR}/unpaired"

RAW_FASTQC="${REPO_ROOT}/results/qc/raw_fastqc"
RAW_MULTIQC="${REPO_ROOT}/results/qc/raw_multiqc"

TRIMMED_FASTQC="${REPO_ROOT}/results/qc/trimmed_fastqc"
TRIMMED_MULTIQC="${REPO_ROOT}/results/qc/trimmed_multiqc"


# ------------------------------------------------------------
# Check required software
# ------------------------------------------------------------

for cmd in trimmomatic fastqc multiqc
do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if ! compgen -G "${RAW_DIR}/*_R1.fastq.gz" > /dev/null; then
    echo "[ERROR] No R1 FASTQ files found in ${RAW_DIR}" >&2
    exit 1
fi


# ------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------

mkdir -p \
    "${TRIMMED_DIR}" \
    "${UNPAIRED_DIR}" \
    "${RAW_FASTQC}" \
    "${RAW_MULTIQC}" \
    "${TRIMMED_FASTQC}" \
    "${TRIMMED_MULTIQC}"


echo "[INFO] Threads: ${THREADS}"
echo "[INFO] Raw reads: ${RAW_DIR}"
echo "[INFO] Trimmed reads: ${TRIMMED_DIR}"


# ------------------------------------------------------------
# 1. Quality control of raw reads
# ------------------------------------------------------------

echo
echo "[INFO] Running FastQC on raw reads"

fastqc \
    "${RAW_DIR}"/*.fastq.gz \
    --threads "${THREADS}" \
    --outdir "${RAW_FASTQC}"


echo
echo "[INFO] Generating MultiQC report for raw reads"

multiqc \
    "${RAW_FASTQC}" \
    --outdir "${RAW_MULTIQC}"


# ------------------------------------------------------------
# 2. Trimming with Trimmomatic
# ------------------------------------------------------------

echo
echo "[INFO] Running Trimmomatic"

for R1 in "${RAW_DIR}"/*_R1.fastq.gz
do

    SAMPLE="$(basename "${R1}" _R1.fastq.gz)"
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "${R2}" ]]; then
        echo "[ERROR] Missing R2 file for ${SAMPLE}" >&2
        exit 1
    fi

    echo "[INFO] Processing ${SAMPLE}"

    trimmomatic PE \
        -threads "${THREADS}" \
        -trimlog "${TRIMMED_DIR}/${SAMPLE}.trimmomatic.log" \
        "${R1}" \
        "${R2}" \
        "${TRIMMED_DIR}/${SAMPLE}_R1.fastq.gz" \
        "${UNPAIRED_DIR}/${SAMPLE}_R1.unpaired.fastq.gz" \
        "${TRIMMED_DIR}/${SAMPLE}_R2.fastq.gz" \
        "${UNPAIRED_DIR}/${SAMPLE}_R2.unpaired.fastq.gz" \
        LEADING:30 \
        TRAILING:30 \
        SLIDINGWINDOW:10:30 \
        MINLEN:50

done


# ------------------------------------------------------------
# 3. Quality control of trimmed paired reads
# ------------------------------------------------------------

echo
echo "[INFO] Running FastQC on trimmed paired reads"

fastqc \
    "${TRIMMED_DIR}"/*_R1.fastq.gz \
    "${TRIMMED_DIR}"/*_R2.fastq.gz \
    --threads "${THREADS}" \
    --outdir "${TRIMMED_FASTQC}"


echo
echo "[INFO] Generating MultiQC report for trimmed reads"

multiqc \
    "${TRIMMED_FASTQC}" \
    --outdir "${TRIMMED_MULTIQC}"


# ------------------------------------------------------------
# Done
# ------------------------------------------------------------

echo
echo "[SUCCESS] Quality control and trimming completed."
echo "[INFO] Trimmed paired reads:"
echo "       ${TRIMMED_DIR}"
echo
echo "[INFO] Raw MultiQC report:"
echo "       ${RAW_MULTIQC}/multiqc_report.html"
echo
echo "[INFO] Trimmed MultiQC report:"
echo "       ${TRIMMED_MULTIQC}/multiqc_report.html"