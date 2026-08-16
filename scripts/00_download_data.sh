#!/usr/bin/env bash
set -euo pipefail

# Download and prepare raw Illumina reads for the H1N2 study.
#
# Data source:
#   NCBI BioProject PRJNA994299
#
# Requirements:
#   - SRA Toolkit (prefetch, fasterq-dump)
#   - gzip
#
# Usage:
#   bash scripts/00_download_data.sh
#
# Optional:
#   THREADS=8 bash scripts/00_download_data.sh
#
# Output:
#   data/fastq/<sample>_R1.fastq.gz
#   data/fastq/<sample>_R2.fastq.gz

THREADS="${THREADS:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
METADATA="${REPO_ROOT}/metadata/samples.tsv"
OUTDIR="${REPO_ROOT}/data/fastq"

if [[ ! -f "${METADATA}" ]]; then
    echo "[ERROR] Metadata file not found: ${METADATA}" >&2
    exit 1
fi

for cmd in prefetch fasterq-dump gzip; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
done

mkdir -p "${OUTDIR}"

echo "[INFO] BioProject: PRJNA994299"
echo "[INFO] Metadata:   ${METADATA}"
echo "[INFO] Output:     ${OUTDIR}"
echo "[INFO] Threads:    ${THREADS}"

# Read the TSV while skipping the header.
tail -n +2 "${METADATA}" | while IFS=$'\t' read -r RUN BIOSAMPLE SAMPLE ORIGINAL_NAME ROLE ISOLATION_SOURCE STRAIN
do
    R1="${OUTDIR}/${SAMPLE}_R1.fastq.gz"
    R2="${OUTDIR}/${SAMPLE}_R2.fastq.gz"

    if [[ -s "${R1}" && -s "${R2}" ]]; then
        echo "[SKIP] ${RUN} -> ${SAMPLE}: FASTQ files already exist."
        continue
    fi

    echo
    echo "[INFO] Processing ${RUN} (${BIOSAMPLE}) -> ${SAMPLE} [${ROLE}]"

    # Cache the SRA run locally.
    prefetch "${RUN}"

    # fasterq-dump writes <RUN>_1.fastq and <RUN>_2.fastq.
    fasterq-dump "${RUN}" \
        --split-files \
        --threads "${THREADS}" \
        --outdir "${OUTDIR}"

    RAW_R1="${OUTDIR}/${RUN}_1.fastq"
    RAW_R2="${OUTDIR}/${RUN}_2.fastq"

    if [[ ! -s "${RAW_R1}" || ! -s "${RAW_R2}" ]]; then
        echo "[ERROR] Expected paired FASTQ files were not created for ${RUN}." >&2
        exit 1
    fi

    mv "${RAW_R1}" "${OUTDIR}/${SAMPLE}_R1.fastq"
    mv "${RAW_R2}" "${OUTDIR}/${SAMPLE}_R2.fastq"

    gzip -f "${OUTDIR}/${SAMPLE}_R1.fastq"
    gzip -f "${OUTDIR}/${SAMPLE}_R2.fastq"

    echo "[DONE] ${RUN} -> ${SAMPLE}"
done

echo
echo "[INFO] Download complete."

EXPECTED=$(($(tail -n +2 "${METADATA}" | wc -l | tr -d ' ') * 2))
FOUND=$(find "${OUTDIR}" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l | tr -d ' ')

echo "[INFO] Expected FASTQ files: ${EXPECTED}"
echo "[INFO] Found FASTQ files:    ${FOUND}"

if [[ "${FOUND}" -ne "${EXPECTED}" ]]; then
    echo "[WARNING] FASTQ count differs from expected value." >&2
    exit 1
fi

echo "[SUCCESS] All paired-end FASTQ files are available in ${OUTDIR}."