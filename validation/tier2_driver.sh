#!/usr/bin/env bash
# Tier 2 driver — disk-managed end-to-end whole-genome trio benchmark.
#
# Sequences download + run + free per sample to stay within a tight disk
# budget (~127 GB free). Each sample takes ~40 min download + ~3 h
# compute + ~30 min hap.py = ~4 h. Sequential 3 samples ≈ ~12 h total.
#
# Disk peak per sample (during compute):
#   BAM (40 GB) + chr1 examples chunk (~50 GB) + reference (3 GB) ≈ 93 GB
#
# Usage:
#   ./validation/tier2_driver.sh                    # full run
#   ./validation/tier2_driver.sh --download-only    # download phase only
#
# Run in background:
#   nohup ./validation/tier2_driver.sh > /tmp/tier2.log 2>&1 &

set -euo pipefail
cd "$(dirname "$0")/.."

DATA="${DV_GIAB_DIR:-/tmp/dv_giab}/full"
GS_BASE="https://storage.googleapis.com/deepvariant/case-study-testdata"
GIAB_FTP="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"

mkdir -p "${DATA}"

dl() {
  local url="$1" out="$2"
  if [ -f "${out}" ] && [ -s "${out}" ]; then
    echo "    ${out} already present, skipping"
    return
  fi
  echo "    Downloading ${out}"
  echo "      URL: ${url}"
  curl -L --fail -C - --progress-bar -o "${out}" "${url}" \
    || { echo "ERROR: download failed for ${url}"; return 1; }
}

download_shared() {
  echo
  echo "=== Stage 1: Reference + truth sets (one-time, small) ==="
  cd "${DATA}"
  # NCBI canonical no_alt_analysis_set reference (matches Google's case-
  # study fixture). Compressed; ~900 MB compressed → ~3.1 GB uncompressed.
  local NCBI_REF="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  dl "${NCBI_REF}" "GRCh38.fa.gz"
  if [ ! -f "GRCh38.fa" ]; then
    echo "    Decompressing GRCh38.fa.gz …"
    gunzip -k "GRCh38.fa.gz"
  fi
  if [ ! -f "GRCh38.fa.fai" ]; then
    echo "    Indexing with samtools faidx …"
    /opt/homebrew/bin/samtools faidx GRCh38.fa
  fi
  for SAMPLE in HG002 HG003 HG004; do
    dl "${GIAB_FTP}/release/AshkenazimTrio/${SAMPLE}_NA*/NISTv4.2.1/GRCh38/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"     "${SAMPLE}.truth.vcf.gz" || true
    dl "${GIAB_FTP}/release/AshkenazimTrio/${SAMPLE}_NA*/NISTv4.2.1/GRCh38/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" "${SAMPLE}.truth.vcf.gz.tbi" || true
    dl "${GIAB_FTP}/release/AshkenazimTrio/${SAMPLE}_NA*/NISTv4.2.1/GRCh38/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "${SAMPLE}.truth.bed" || true
  done
  cd - > /dev/null

  # Truth sets fall back to /tmp/dv_giab/data/ for backward compatibility.
  for s in HG002 HG003 HG004; do
    if [ ! -f "${DATA}/${s}.truth.vcf.gz" ] && [ -f "/tmp/dv_giab/data/${s}.truth.vcf.gz" ]; then
      ln -sf "/tmp/dv_giab/data/${s}.truth.vcf.gz" "${DATA}/${s}.truth.vcf.gz"
      ln -sf "/tmp/dv_giab/data/${s}.truth.vcf.gz.tbi" "${DATA}/${s}.truth.vcf.gz.tbi"
      ln -sf "/tmp/dv_giab/data/${s}.truth.bed" "${DATA}/${s}.truth.bed"
    fi
  done
  # HG002 truth.vcf.gz (plain name, legacy) → HG002.truth.vcf.gz.
  if [ ! -f "${DATA}/HG002.truth.vcf.gz" ] && [ -f "/tmp/dv_giab/data/truth.vcf.gz" ]; then
    ln -sf "/tmp/dv_giab/data/truth.vcf.gz" "${DATA}/HG002.truth.vcf.gz"
    ln -sf "/tmp/dv_giab/data/truth.vcf.gz.tbi" "${DATA}/HG002.truth.vcf.gz.tbi"
    ln -sf "/tmp/dv_giab/data/truth.bed" "${DATA}/HG002.truth.bed"
  fi
}

process_one_sample() {
  local sample="$1"
  echo
  echo "============================================================"
  echo "=== ${sample} : download → run → hap.py → free BAM"
  echo "============================================================"
  date -u +"    Started: %Y-%m-%dT%H:%M:%SZ"

  # Disk pre-check.
  local free_gb
  free_gb=$(df -g "${DATA}" | tail -1 | awk '{print $4}')
  echo "    Disk free at start: ${free_gb} GB"
  if [ "${free_gb}" -lt 50 ]; then
    echo "ERROR: only ${free_gb} GB free — need ≥50 GB headroom for chunk processing" >&2
    return 1
  fi

  # Stage 1: download BAM (~40 GB, ~30-60 min).
  echo
  echo "--- ${sample}: download BAM ---"
  dl "${GS_BASE}/${sample}.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"     "${DATA}/${sample}.bam"
  dl "${GS_BASE}/${sample}.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai" "${DATA}/${sample}.bam.bai"

  # Stage 2: chunked WG run + hap.py.
  echo
  echo "--- ${sample}: chunked WG run + hap.py ---"
  DV_GIAB_DIR="${DV_GIAB_DIR:-/tmp/dv_giab}" \
    ./validation/run_giab_wg_chunked.sh "${sample}"
}

# ── Main ─────────────────────────────────────────────────────────────────
echo "Tier 2 driver — whole-genome trio benchmark"
echo "Working dir : ${DATA}"
date -u +"Started     : %Y-%m-%dT%H:%M:%SZ"

download_shared

if [ "${1:-}" = "--download-only" ]; then
  echo "==> --download-only requested; reference + truth sets fetched"
  exit 0
fi

# Process each sample sequentially.
process_one_sample HG002
process_one_sample HG003
process_one_sample HG004

# Final summary already produced by run_giab_wg_chunked.sh.
echo
date -u +"Finished    : %Y-%m-%dT%H:%M:%SZ"
echo
echo "All done. Inspect:"
echo "  validation/output/wg_trio_summary.tsv"
echo "  validation/output/{HG002,HG003,HG004}_wg/happy.summary.csv"
