#!/usr/bin/env bash
# Download GIAB genome stratifications v3.6 for stratified hap.py runs.
#
# Stratification BED files split GIAB confidence regions into context
# categories (LowComplexity, SegmentalDuplications, MHC, GC bands, etc.).
# Used by hap.py via the `--stratification <tsv>` flag to produce
# per-context F1 breakdowns. Tells us WHERE our remaining FN/FP sit —
# lowcomplexity vs segdup vs MHC — and guides which model improvements
# would give the biggest F1 lift.
#
# Source: GIAB FTP, v3.6 GRCh38 (Apr 2025).
# Size:   ~1.4 GB compressed; ~5 GB extracted.
#
# Usage:
#   ./validation/download_giab_strats.sh [target_dir]
#     target_dir defaults to ${DV_GIAB_DIR}/strats or /tmp/dv_giab/strats.
#
# After download, run a stratified hap.py:
#   docker run --rm --platform linux/amd64 \
#     -v /tmp/dv_giab/data:/data:ro \
#     -v /tmp/dv_giab/strats:/strats:ro \
#     -v /path/to/giab_test:/work \
#     jmcdani20/hap.py:v0.3.12 \
#     /opt/hap.py/bin/hap.py /data/truth.vcf.gz /work/our.vcf.gz \
#       -f /data/truth.bed -r /data/GRCh38.fa \
#       -o /work/happy_strat \
#       --location chr20 \
#       --stratification /strats/v3.6-stratifications-GRCh38.tsv

set -euo pipefail

TARGET="${1:-${DV_GIAB_DIR:-/tmp/dv_giab}/strats}"
URL="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.6/genome-stratifications-GRCh38@all.tar.gz"

mkdir -p "${TARGET}"
cd "${TARGET}"

if [ -f "v3.6-stratifications-GRCh38.tsv" ]; then
  echo "==> Stratifications already present at ${TARGET} — skipping download"
  exit 0
fi

ARCHIVE="genome-stratifications-GRCh38@all.tar.gz"
if [ ! -f "${ARCHIVE}" ]; then
  echo "==> Downloading GIAB stratifications v3.6 GRCh38 (~1.4 GB) ..."
  echo "    URL: ${URL}"
  echo "    Target: ${TARGET}/${ARCHIVE}"
  curl -fL --retry 3 --connect-timeout 15 --progress-bar -o "${ARCHIVE}.partial" "${URL}"
  mv "${ARCHIVE}.partial" "${ARCHIVE}"
fi

echo "==> Extracting ..."
tar xzf "${ARCHIVE}"

# The TSV index after extraction is typically at v3.6-stratifications-GRCh38.tsv
if [ ! -f "v3.6-stratifications-GRCh38.tsv" ]; then
  # Try common alternative names.
  TSV=$(find . -name "*stratifications*GRCh38*.tsv" -type f | head -1)
  if [ -n "${TSV}" ]; then
    ln -sf "${TSV}" v3.6-stratifications-GRCh38.tsv
  else
    echo "WARNING: could not locate the .tsv index. Inspect ${TARGET}/" >&2
  fi
fi

echo
echo "==> Done. Stratifications at: ${TARGET}/"
echo "    TSV index: ${TARGET}/v3.6-stratifications-GRCh38.tsv"
echo
echo "    Total size:"
du -sh "${TARGET}/" 2>/dev/null
