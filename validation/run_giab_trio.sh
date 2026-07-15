#!/usr/bin/env bash
# Phase 9 / Step 5b — whole-genome F1 trio validation orchestrator.
#
# Prereq: ./validation/download_giab_full_genome.sh has been run (~120 GB
# of data at ${DV_GIAB_DIR}/full).
#
# Runs deepvariant + hap.py for HG002, HG003, HG004 in sequence
# (~3 hours per sample on M4 Max 14-thread, ~9 hours total). Captures
# F1 metrics into validation/output/<sample>_full_genome/happy.summary.csv.
#
# Phase 4 spec gates (per sample, on whole-genome GIAB high-confidence
# regions):
#   SNP F1   ≥ Linux upstream F1 − 0.05 %
#   INDEL F1 ≥ Linux upstream F1 − 0.10 %
#
# Usage:
#   ./validation/run_giab_trio.sh
#
# Each sample is independent — failure on HG003 doesn't block HG004.

set -euo pipefail
cd "$(dirname "$0")/.."

DATA="${DV_GIAB_DIR:-/tmp/dv_giab}/full"

if [ ! -f "${DATA}/HG002.bam" ] || [ ! -f "${DATA}/GRCh38.fa" ]; then
  echo "ERROR: missing data at ${DATA}" >&2
  echo "  Run ./validation/download_giab_full_genome.sh first." >&2
  exit 1
fi

run_sample() {
  local sample="$1" truth_vcf="$2" truth_bed="$3"
  local out="validation/output/${sample}_full_genome"
  mkdir -p "${out}"

  echo
  echo "============================================================"
  echo "==> ${sample} whole-genome (~3 h on M4 Max 14-thread)"
  echo "============================================================"
  echo "    BAM:   ${DATA}/${sample}.bam"
  echo "    Truth: ${truth_vcf}"
  echo "    BED:   ${truth_bed}"
  echo

  # Stage 1: deepvariant native run.
  if [ ! -f "${out}/our.vcf.gz" ]; then
    DV_GIAB_DIR="${DATA}" DV_VALIDATION_OUT="${out}" \
      ./validation/run_giab.sh "${sample}" ""    # empty region = whole genome
  else
    echo "==> ${out}/our.vcf.gz exists, skipping deepvariant run"
  fi

  # Stage 2: hap.py vs sample's truth set.
  if [ ! -f "${out}/happy.summary.csv" ]; then
    echo
    echo "==> hap.py vs ${sample} GIAB v4.2.1 truth (whole-genome)"
    docker run --rm --platform linux/amd64 \
      -v "${DATA}:/data:ro" \
      -v "$(realpath "${out}"):/work" \
      jmcdani20/hap.py:v0.3.12 \
      /opt/hap.py/bin/hap.py \
        "/data/$(basename "${truth_vcf}")" \
        /work/our.vcf.gz \
        -f "/data/$(basename "${truth_bed}")" \
        -r /data/GRCh38.fa \
        -o /work/happy
  else
    echo "==> ${out}/happy.summary.csv exists, skipping hap.py"
  fi

  echo
  echo "==> ${sample} F1 (PASS rows):"
  awk -F, 'NR==1 || /,PASS,/ {print}' "${out}/happy.summary.csv"
}

# HG002: truth already at /tmp/dv_giab/data/truth.vcf.gz (legacy location).
run_sample HG002 "${DATA}/HG002.truth.vcf.gz" "${DATA}/HG002.truth.bed" \
  || run_sample HG002 "/tmp/dv_giab/data/truth.vcf.gz" "/tmp/dv_giab/data/truth.bed"

# HG003 + HG004: downloaded by download_giab_full_genome.sh.
run_sample HG003 "${DATA}/HG003.truth.vcf.gz" "${DATA}/HG003.truth.bed"
run_sample HG004 "${DATA}/HG004.truth.vcf.gz" "${DATA}/HG004.truth.bed"

echo
echo "============================================================"
echo "==> Trio whole-genome F1 validation complete."
echo "    Results at validation/output/{HG002,HG003,HG004}_full_genome/"
echo "============================================================"
