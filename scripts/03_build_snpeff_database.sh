#!/usr/bin/env bash
set -euo pipefail

# Build a local custom SnpEff database for the H1N2 reference.
#
# This script reconstructs the procedure used in the original analysis,
# but keeps the SnpEff database entirely inside the repository.
#
# Requirements:
#   snpEff
#   gffread
#
# Input:
#   reference/H1N2_reference.fasta
#   reference/H1N2_reference.gff3
#
# Output:
#   resources/snpeff/data/H1N2/
#
# Usage:
#   bash scripts/03_build_snpeff_database.sh


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DB_NAME="H1N2"

REFERENCE="${REPO_ROOT}/reference/H1N2_reference.fasta"
GFF3="${REPO_ROOT}/reference/H1N2_reference.gff3"

SNPEFF_ROOT="${REPO_ROOT}/resources/snpeff"
SNPEFF_DATA="${SNPEFF_ROOT}/data"
DB_DIR="${SNPEFF_DATA}/${DB_NAME}"

CONFIG="${SNPEFF_ROOT}/snpEff.config"
GTF="${DB_DIR}/genes.gtf"
SEQUENCES="${DB_DIR}/sequences.fa"


# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

for cmd in snpEff gffread
do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

for FILE in "${REFERENCE}" "${GFF3}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "[ERROR] Missing or empty input: ${FILE}" >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# Create local SnpEff structure
# ------------------------------------------------------------

mkdir -p "${DB_DIR}"


# ------------------------------------------------------------
# Prepare GFF3 with segment names used in the workflow
# ------------------------------------------------------------

echo "[INFO] Normalizing segment identifiers in GFF3"

NORMALIZED_GFF3="${DB_DIR}/H1N2.normalized.gff3"

awk 'BEGIN{FS=OFS="\t"}
{
    if ($0 ~ /^#/) {
        print;
        next;
    }

    if ($1=="HF674912.1") $1="HA";
    else if ($1=="HF674913.1") $1="NA";
    else if ($1=="HF674914.1") $1="M";
    else if ($1=="HF674915.1") $1="NP";
    else if ($1=="HF674916.1") $1="NS";
    else if ($1=="HF674917.1") $1="PA";
    else if ($1=="HF674918.1") $1="PB1";
    else if ($1=="HF674919.1") $1="PB2";

    print;
}' "${GFF3}" > "${NORMALIZED_GFF3}"


# ------------------------------------------------------------
# Convert GFF3 to GTF
# ------------------------------------------------------------

echo "[INFO] Converting GFF3 to GTF"

gffread \
    "${NORMALIZED_GFF3}" \
    -T \
    -o "${GTF}"


# ------------------------------------------------------------
# Copy reference FASTA
# ------------------------------------------------------------

echo "[INFO] Preparing reference FASTA"

cp "${REFERENCE}" "${SEQUENCES}"


# ------------------------------------------------------------
# Create local SnpEff config
# ------------------------------------------------------------

echo "[INFO] Creating local SnpEff configuration"

cat > "${CONFIG}" <<EOF
data.dir = ${SNPEFF_DATA}

${DB_NAME}.genome : Swine influenza H1N2
EOF


# ------------------------------------------------------------
# Build database
# ------------------------------------------------------------

echo "[INFO] Building SnpEff database: ${DB_NAME}"

snpEff build \
    -gtf22 \
    -v \
    -noCheckCds \
    -noCheckProtein \
    -c "${CONFIG}" \
    "${DB_NAME}"


# ------------------------------------------------------------
# Validate database
# ------------------------------------------------------------

echo "[INFO] Validating SnpEff database"

SNPEFF_BIN="${DB_DIR}/snpEffectPredictor.bin"

if [[ ! -s "${SNPEFF_BIN}" ]]; then
    echo "[ERROR] SnpEff database binary was not created: ${SNPEFF_BIN}" >&2
    exit 1
fi


echo
echo "[SUCCESS] SnpEff database built successfully."
echo "[INFO] Database: ${DB_NAME}"
echo "[INFO] Config:   ${CONFIG}"
echo "[INFO] Data dir: ${DB_DIR}"