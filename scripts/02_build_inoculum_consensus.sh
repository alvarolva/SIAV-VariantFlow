#!/usr/bin/env bash
set -euo pipefail

# Build the inoculum consensus sequence
# Reconstructed from the historical EVOLSI workflow.
#
# Input:
#   reference/H1N2_reference.fasta
#   data/trimmed/H1N2_R1.fastq.gz
#   data/trimmed/H1N2_R2.fastq.gz
#
# Output:
#   results/consensus/H1N2_inoculum_consensus.fasta
#
# Requirements:
#   bwa
#   samtools
#   gatk4
#   bcftools
#
# Usage:
#   bash scripts/02_build_inoculum_consensus.sh
#
# Optional:
#   THREADS=8 bash scripts/02_build_inoculum_consensus.sh


THREADS="${THREADS:-4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REFERENCE="${REPO_ROOT}/reference/H1N2_reference.fasta"

R1="${REPO_ROOT}/data/trimmed/H1N2_R1.fastq.gz"
R2="${REPO_ROOT}/data/trimmed/H1N2_R2.fastq.gz"

OUTDIR="${REPO_ROOT}/results/consensus"
SAMPLE="H1N2"

SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
CLEAN_BAM="${OUTDIR}/${SAMPLE}.clean.sorted.bam"
DEDUP_BAM="${OUTDIR}/${SAMPLE}.dedup.bam"

VCF="${OUTDIR}/${SAMPLE}.bcftools.vcf.gz"

CONSENSUS="${OUTDIR}/H1N2_inoculum_consensus.fasta"


# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

for cmd in bwa samtools gatk bcftools
do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

for FILE in "${REFERENCE}" "${R1}" "${R2}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "[ERROR] Missing or empty input: ${FILE}" >&2
        exit 1
    fi
done


mkdir -p "${OUTDIR}"

echo "[INFO] Threads:   ${THREADS}"
echo "[INFO] Reference: ${REFERENCE}"
echo "[INFO] R1:        ${R1}"
echo "[INFO] R2:        ${R2}"


# ------------------------------------------------------------
# 1. Index reference
# ------------------------------------------------------------

echo
echo "[INFO] Indexing reference"

if [[ ! -f "${REFERENCE}.bwt" ]]; then
    bwa index "${REFERENCE}"
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    samtools faidx "${REFERENCE}"
fi


# ------------------------------------------------------------
# 2. Align inoculum reads
# ------------------------------------------------------------

echo
echo "[INFO] Aligning inoculum reads with BWA-MEM"

READ_GROUP="@RG\tID:H1N2\tSM:H1N2\tLB:nextera\tPL:ILLUMINA"

bwa mem \
    -M \
    -t "${THREADS}" \
    -R "${READ_GROUP}" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${SORTED_BAM}"

samtools index "${SORTED_BAM}"


# ------------------------------------------------------------
# 3. Remove unmapped reads and MAPQ <30
# ------------------------------------------------------------

echo
echo "[INFO] Filtering unmapped reads and MAPQ <30"

samtools view \
    -h \
    -F 4 \
    -b \
    -q 30 \
    "${SORTED_BAM}" \
    > "${CLEAN_BAM}"

samtools index "${CLEAN_BAM}"


# ------------------------------------------------------------
# 4. Remove sequencing duplicates
# ------------------------------------------------------------

echo
echo "[INFO] Removing sequencing duplicates"

gatk MarkDuplicatesSpark \
    -I "${CLEAN_BAM}" \
    -O "${DEDUP_BAM}" \
    --remove-sequencing-duplicates \
    -M "${OUTDIR}/${SAMPLE}.dedup.metrics.txt"

samtools index "${DEDUP_BAM}"


# ------------------------------------------------------------
# 5. Alignment QC
# ------------------------------------------------------------

echo
echo "[INFO] Generating alignment statistics"

samtools flagstat \
    "${DEDUP_BAM}" \
    > "${OUTDIR}/${SAMPLE}.flagstat.txt"

samtools idxstats \
    "${DEDUP_BAM}" \
    > "${OUTDIR}/${SAMPLE}.idxstats.txt"


# ------------------------------------------------------------
# 6. Call consensus-level variants with bcftools
# ------------------------------------------------------------

echo
echo "[INFO] Calling variants for consensus construction"

bcftools mpileup \
    -Ou \
    -f "${REFERENCE}" \
    "${DEDUP_BAM}" \
    | bcftools call \
        -mv \
        -Oz \
        -o "${VCF}"

bcftools index -f "${VCF}"


# ------------------------------------------------------------
# 7. Build inoculum consensus
# ------------------------------------------------------------

echo
echo "[INFO] Building inoculum consensus"

bcftools consensus \
    -f "${REFERENCE}" \
    "${VCF}" \
    > "${CONSENSUS}"


# ------------------------------------------------------------
# 8. Validate consensus
# ------------------------------------------------------------

echo
echo "[INFO] Validating consensus"

samtools faidx "${CONSENSUS}"

SEGMENT_COUNT=$(grep -c '^>' "${CONSENSUS}")

if [[ "${SEGMENT_COUNT}" -ne 8 ]]; then
    echo "[ERROR] Expected 8 segments, found ${SEGMENT_COUNT}" >&2
    exit 1
fi

EXPECTED_SEGMENTS="HA NA M NP NS PA PB1 PB2"

for SEGMENT in ${EXPECTED_SEGMENTS}
do
    if ! grep -Fxq ">${SEGMENT}" "${CONSENSUS}"; then
        echo "[ERROR] Missing expected segment: ${SEGMENT}" >&2
        exit 1
    fi
done


echo
echo "[SUCCESS] Inoculum consensus generated successfully."
echo
echo "[INFO] Consensus:"
echo "       ${CONSENSUS}"
echo
echo "[INFO] Segment lengths:"
cut -f1,2 "${CONSENSUS}.fai"