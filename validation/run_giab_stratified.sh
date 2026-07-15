#!/usr/bin/env bash
# Tier 3 — stratified F1 breakdown via GIAB stratifications v3.6.
#
# Re-runs hap.py against an already-evaluated VCF with a stratification
# TSV. Output: per-context F1 (LowComplexity, SegmentalDuplications,
# MHC, GC bands, Functional, etc.) — tells us WHERE our remaining
# FN/FP sit, guiding which model improvements would yield the biggest
# F1 lift.
#
# Prerequisites:
#   - validation/output/<sample>_<region>/our.vcf.gz already exists
#   - /tmp/dv_giab/strats/ populated by validation/download_giab_strats.sh
#
# Usage:
#   ./validation/run_giab_stratified.sh HG002 chr20
#   ./validation/run_giab_stratified.sh HG002 wg
#   ./validation/run_giab_stratified.sh HG003 chr20
#
# Wall-time: ~10-15 min per sample (hap.py runs longer with stratification)

set -euo pipefail
cd "$(dirname "$0")/.."

SAMPLE="${1:?usage: $0 <sample> <region>}"
REGION="${2:?usage: $0 <sample> <region>  region in {chr20, wg}}"

DATA="${DV_GIAB_DIR:-/tmp/dv_giab/data}"
STRATS="${DV_STRATS_DIR:-/tmp/dv_giab/strats/GRCh38@all}"

if [ ! -f "${STRATS}/GRCh38-all-stratifications.tsv" ]; then
  echo "ERROR: stratifications missing at ${STRATS}" >&2
  echo "       Run ./validation/download_giab_strats.sh first." >&2
  exit 1
fi

# Determine OUT path + truth path based on region.
case "${REGION}" in
  chr20) OUT="validation/output/${SAMPLE}_chr20" ;;
  wg)    OUT="validation/output/${SAMPLE}_wg" ;;
  *)     echo "ERROR: region must be 'chr20' or 'wg' (got ${REGION})"; exit 1 ;;
esac

# Truth set per sample. HG002's legacy filename is `truth.vcf.gz`,
# others use HG00x.truth.vcf.gz. The Tier 2 driver normalises this.
case "${SAMPLE}" in
  HG002)
    TRUTH_VCF="${DATA}/truth.vcf.gz"
    TRUTH_BED="${DATA}/truth.bed"
    [ -f "${DATA}/HG002.truth.vcf.gz" ] && TRUTH_VCF="${DATA}/HG002.truth.vcf.gz"
    [ -f "${DATA}/HG002.truth.bed" ] && TRUTH_BED="${DATA}/HG002.truth.bed"
    ;;
  HG003|HG004)
    TRUTH_VCF="${DATA}/${SAMPLE}.truth.vcf.gz"
    TRUTH_BED="${DATA}/${SAMPLE}.truth.bed"
    ;;
  *)
    echo "ERROR: unknown sample ${SAMPLE}"; exit 1 ;;
esac

if [ ! -f "${OUT}/our.vcf.gz" ]; then
  echo "ERROR: ${OUT}/our.vcf.gz missing — run the per-region pipeline first." >&2
  exit 1
fi

echo "=========================================================="
echo "==> Stratified F1: ${SAMPLE} ${REGION}"
echo "    VCF:   ${OUT}/our.vcf.gz"
echo "    Truth: ${TRUTH_VCF}"
echo "    BED:   ${TRUTH_BED}"
echo "    Strats: ${STRATS}/GRCh38-all-stratifications.tsv"
echo "=========================================================="

LOC_ARGS=()
[ "${REGION}" = "chr20" ] && LOC_ARGS+=(--location chr20)

docker run --rm \
  -v "${DATA}:/data:ro" \
  -v "${STRATS}:/strats:ro" \
  -v "$(realpath "${OUT}"):/work" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    "/data/$(basename "${TRUTH_VCF}")" \
    /work/our.vcf.gz \
    -f "/data/$(basename "${TRUTH_BED}")" \
    -r /data/GRCh38.fa \
    -o /work/happy_strat \
    --stratification /strats/GRCh38-all-stratifications.tsv \
    "${LOC_ARGS[@]}"

echo
echo "==> Stratified F1 saved at ${OUT}/happy_strat.summary.csv"
echo "    Inspect by-context breakdown:"
echo "      column -t -s, ${OUT}/happy_strat.extended.csv | less -S"
echo
echo "==> Top-10 contexts by F1 deficit (PASS, SNP):"
awk -F, 'NR>1 && $2=="PASS" && $1=="SNP" && $5 != "*" && $14 != "" && $14+0 < 1.0 {
           printf "%-50s F1=%-8s n=%s\n", $5, $14, $3
         }' "${OUT}/happy_strat.extended.csv" 2>/dev/null \
  | sort -k2,2 \
  | head -10 || true
