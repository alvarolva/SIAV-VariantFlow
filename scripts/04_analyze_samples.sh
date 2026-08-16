#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 04_analyze_samples.sh
#
# Mapping, coverage and low-frequency variant analysis.
# Reconstructed from the published H1N2 influenza workflow.
#
# INPUT:
#   data/trimmed/<sample>_R1.fastq.gz
#   data/trimmed/<sample>_R2.fastq.gz
#
# REFERENCE:
#   results/consensus/H1N2_inoculum_consensus.fasta
#
# USAGE:
#   Single sample:
#       bash scripts/04_analyze_samples.sh 01_BALF
#
#   All samples:
#       bash scripts/04_analyze_samples.sh --all
#
# Optional:
#       THREADS=8 bash scripts/04_analyze_samples.sh --all
#
# Published SNV criteria:
#   minimum depth             = 100 reads
#   minimum alternative count = 10 reads
#
# LoFreq performs its own statistical variant-quality filtering.
# ============================================================

THREADS="${THREADS:-4}"
MIN_DEPTH="${MIN_DEPTH:-100}"
MIN_ALT_COUNT="${MIN_ALT_COUNT:-10}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

FASTQ_DIR="${REPO_ROOT}/data/trimmed"

REFERENCE="${REPO_ROOT}/results/consensus/H1N2_inoculum_consensus.fasta"

RESULTS_DIR="${REPO_ROOT}/results/variants"
SUMMARY_DIR="${REPO_ROOT}/results/summary"

SNPEFF_CONFIG="${REPO_ROOT}/resources/snpeff/snpEff.config"
SNPEFF_DB="H1N2"

BOWTIE_INDEX="${REPO_ROOT}/results/reference_index/H1N2_inoculum"


# ============================================================
# Requirements
# ============================================================

for cmd in \
    bowtie2 \
    bowtie2-build \
    samtools \
    bcftools \
    gatk \
    lofreq \
    snpEff
do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
done


# ============================================================
# Input checks
# ============================================================

if [[ ! -s "${REFERENCE}" ]]; then
    echo "[ERROR] Consensus reference not found: ${REFERENCE}" >&2
    exit 1
fi

if [[ ! -s "${SNPEFF_CONFIG}" ]]; then
    echo "[ERROR] SnpEff config not found: ${SNPEFF_CONFIG}" >&2
    echo "[INFO] Run scripts/03_build_snpeff_database.sh first." >&2
    exit 1
fi

mkdir -p \
    "${RESULTS_DIR}" \
    "${SUMMARY_DIR}" \
    "$(dirname "${BOWTIE_INDEX}")"


# ============================================================
# Reference preparation
# ============================================================

echo "[INFO] Preparing reference"

if [[ ! -f "${REFERENCE}.fai" ]]; then
    samtools faidx "${REFERENCE}"
fi

DICT="${REFERENCE%.*}.dict"

if [[ ! -f "${DICT}" ]]; then
    gatk CreateSequenceDictionary \
        -R "${REFERENCE}" \
        -O "${DICT}"
fi


if [[ ! -f "${BOWTIE_INDEX}.1.bt2" && \
      ! -f "${BOWTIE_INDEX}.1.bt2l" ]]; then

    echo "[INFO] Building Bowtie2 index"

    bowtie2-build \
        "${REFERENCE}" \
        "${BOWTIE_INDEX}"
fi


# ============================================================
# Function: analyse one sample
# ============================================================

analyse_sample() {

    SAMPLE="$1"

    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE}"

    mkdir -p "${SAMPLE_DIR}"


    # --------------------------------------------------------
    # Check paired reads
    # --------------------------------------------------------

    if [[ ! -s "${R1}" ]]; then
        echo "[ERROR] Missing R1 for ${SAMPLE}: ${R1}" >&2
        return 1
    fi

    if [[ ! -s "${R2}" ]]; then
        echo "[ERROR] Missing R2 for ${SAMPLE}: ${R2}" >&2
        return 1
    fi


    echo
    echo "============================================================"
    echo "[INFO] Sample: ${SAMPLE}"
    echo "============================================================"


    # --------------------------------------------------------
    # 1. Bowtie2 mapping
    # --------------------------------------------------------

    echo "[INFO] Mapping reads"

    SORTED_BAM="${SAMPLE_DIR}/${SAMPLE}.sorted.bam"

    bowtie2 \
        -x "${BOWTIE_INDEX}" \
        -1 "${R1}" \
        -2 "${R2}" \
        -p "${THREADS}" \
        --very-sensitive-local \
        -a \
        --rg-id "${SAMPLE}" \
        --rg "SM:${SAMPLE}" \
        --rg "LB:nextera" \
        --rg "PL:ILLUMINA" \
        | samtools sort \
            -@ "${THREADS}" \
            -o "${SORTED_BAM}"

    samtools index "${SORTED_BAM}"


    # --------------------------------------------------------
    # 2. Remove unmapped and MAPQ <30
    # --------------------------------------------------------

    echo "[INFO] Filtering unmapped reads and MAPQ <30"

    CLEAN_BAM="${SAMPLE_DIR}/${SAMPLE}.clean.sorted.bam"

    samtools view \
        -h \
        -F 4 \
        -b \
        -q 30 \
        "${SORTED_BAM}" \
        > "${CLEAN_BAM}"

    samtools index "${CLEAN_BAM}"


    # --------------------------------------------------------
    # 3. Generate candidate sites for BQSR
    #
    # Reconstructs the strategy used in the original workflow:
    # candidate variation in each sample was used as known-sites.
    # --------------------------------------------------------

    echo "[INFO] Generating candidate sites for recalibration"

    KNOWN_VCF="${SAMPLE_DIR}/${SAMPLE}.known_sites.vcf.gz"

    bcftools mpileup \
        -Ou \
        -f "${REFERENCE}" \
        "${CLEAN_BAM}" \
        | bcftools call \
            -mv \
            -Oz \
            -o "${KNOWN_VCF}"

    bcftools index -f "${KNOWN_VCF}"

    gatk IndexFeatureFile \
    -I "${KNOWN_VCF}"

    # --------------------------------------------------------
    # 4. Remove sequencing duplicates
    # --------------------------------------------------------

    echo "[INFO] Removing sequencing duplicates"

    DEDUP_BAM="${SAMPLE_DIR}/${SAMPLE}.dedup.bam"

    gatk MarkDuplicatesSpark \
        -I "${CLEAN_BAM}" \
        -O "${DEDUP_BAM}" \
        --remove-sequencing-duplicates \
        -M "${SAMPLE_DIR}/${SAMPLE}.dedup.metrics.txt"

    samtools index "${DEDUP_BAM}"


    # --------------------------------------------------------
    # 5. Base quality recalibration
    # --------------------------------------------------------

    echo "[INFO] Base quality recalibration"

    RECAL_TABLE="${SAMPLE_DIR}/${SAMPLE}.recal_data.table"

    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${DEDUP_BAM}" \
        --known-sites "${KNOWN_VCF}" \
        -O "${RECAL_TABLE}"


    # --------------------------------------------------------
    # 6. Coverage
    #
    # Output per sample:
    # segment    position    coverage    sample
    # --------------------------------------------------------

    echo "[INFO] Calculating coverage"

    DEPTH="${SAMPLE_DIR}/${SAMPLE}.depth.tsv"

    {
        printf "segment\tposition\tcoverage\tsample\n"

        samtools depth \
            -a \
            "${DEDUP_BAM}" \
            | awk -v sample="${SAMPLE}" \
                'BEGIN{OFS="\t"} {print $1,$2,$3,sample}'

    } > "${DEPTH}"


    # --------------------------------------------------------
    # Mapping statistics
    # --------------------------------------------------------

    samtools flagstat \
        "${DEDUP_BAM}" \
        > "${SAMPLE_DIR}/${SAMPLE}.flagstat.txt"


    {
        printf "segment\tmapped_reads\tsample\n"

        samtools idxstats "${DEDUP_BAM}" \
            | awk -v sample="${SAMPLE}" \
                'BEGIN{OFS="\t"} $1!="*" {print $1,$3,sample}'

    } > "${SAMPLE_DIR}/${SAMPLE}.reads_per_segment.tsv"


    # --------------------------------------------------------
    # 7. LoFreq variant calling
    # --------------------------------------------------------

    echo "[INFO] Calling variants with LoFreq"

    RAW_VCF="${SAMPLE_DIR}/${SAMPLE}.lofreq.vcf"

    lofreq call-parallel \
        --pp-threads "${THREADS}" \
        -f "${REFERENCE}" \
        -o "${RAW_VCF}" \
        "${DEDUP_BAM}"


    # --------------------------------------------------------
    # 8. Published filtering criteria
    #
    # DP >= 100
    # alternative reads >= 10
    #
    # LoFreq AF = alternative allele frequency
    # alt count ≈ DP * AF
    # --------------------------------------------------------

    echo "[INFO] Applying published variant filters"

    FILTERED_VCF="${SAMPLE_DIR}/${SAMPLE}.filtered.vcf"

    bcftools filter \
        -i "DP>=${MIN_DEPTH} && DP*AF>=${MIN_ALT_COUNT}" \
        "${RAW_VCF}" \
        > "${FILTERED_VCF}"


    # --------------------------------------------------------
    # 9. SnpEff annotation
    # --------------------------------------------------------

    echo "[INFO] Annotating variants with SnpEff"

    SNPEFF_VCF="${SAMPLE_DIR}/${SAMPLE}.snpeff.vcf"

    conda run -n snpeff_env snpEff \
    -c "${SNPEFF_CONFIG}" \
    -noStats \
    "${SNPEFF_DB}" \
    "${FILTERED_VCF}" \
    > "${SNPEFF_VCF}"


    # --------------------------------------------------------
    # 10. Per-sample table
    #
    # SnpEff ANN annotations are expanded so that each
    # annotation occupies one row.
    # --------------------------------------------------------

    TABLE="${SAMPLE_DIR}/${SAMPLE}.variants.tsv"

    printf \
    "segment\tposition\tREF\tALT\tQUAL\tsample\tDP\tAF\tANN\n" \
    > "${TABLE}"


    bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AF\t%INFO/ANN\n' \
    "${SNPEFF_VCF}" \
    | awk -v sample="${SAMPLE}" \
        'BEGIN{FS=OFS="\t"}
         {
             print $1,$2,$3,$4,$5,sample,$6,$7,$8
         }' \
    >> "${TABLE}"


    echo "[SUCCESS] ${SAMPLE}"
}


# ============================================================
# Decide whether to analyse one sample or all
# ============================================================

if [[ $# -ne 1 ]]; then

    echo "Usage:"
    echo "  bash scripts/04_analyze_samples.sh SAMPLE"
    echo "  bash scripts/04_analyze_samples.sh --all"

    exit 1
fi


if [[ "$1" == "--all" ]]; then

    echo "[INFO] Detecting samples in ${FASTQ_DIR}"

    SAMPLES=()

    for R1 in "${FASTQ_DIR}"/*_R1.fastq.gz
    do
        SAMPLE="$(basename "${R1}" _R1.fastq.gz)"
        SAMPLES+=("${SAMPLE}")
    done

    echo "[INFO] Samples detected: ${#SAMPLES[@]}"

else

    SAMPLES=("$1")

fi


# ============================================================
# Run
# ============================================================

for SAMPLE in "${SAMPLES[@]}"
do
    analyse_sample "${SAMPLE}"
done


# ============================================================
# Combined outputs
# ============================================================

echo
echo "[INFO] Building combined coverage table"

ALL_DEPTH="${SUMMARY_DIR}/all_depth.tsv"

printf "segment\tposition\tcoverage\tsample\n" \
    > "${ALL_DEPTH}"

for SAMPLE in "${SAMPLES[@]}"
do

    tail -n +2 \
        "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}.depth.tsv" \
        >> "${ALL_DEPTH}"

done


echo "[INFO] Building combined variant table"

ALL_VARIANTS="${SUMMARY_DIR}/all_variants.tsv"

printf \
    "segment\tposition\tREF\tALT\tQUAL\tsample\tDP\tAF\tANN\n" \
    > "${ALL_VARIANTS}"

for SAMPLE in "${SAMPLES[@]}"
do

    tail -n +2 \
        "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}.variants.tsv" \
        >> "${ALL_VARIANTS}"

done


echo
echo "============================================================"
echo "[SUCCESS] Analysis completed"
echo "============================================================"
echo
echo "[INFO] Combined coverage:"
echo "       ${ALL_DEPTH}"
echo
echo "[INFO] Combined variants:"
echo "       ${ALL_VARIANTS}"
echo