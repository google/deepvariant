#!/usr/bin/env bash
# Auto-runs the upstream `google/deepvariant:1.10.0` Docker on HG002 WG
# with the SAME configuration as our native run (--num_shards=14), so we
# can do an apples-to-apples bit-level comparison.
#
# Triggered automatically by run_hg002_wg_compare.sh once our native
# HG002 WG completes. Wall-time on this M4 Max ≈ 22 h (Rosetta 2 amd64
# emulation, no GPU passthrough). Run overnight.
#
# Outputs:
#   /tmp/dv_hg002_wg_docker/output.vcf.gz
#   /tmp/dv_hg002_wg_docker/output.g.vcf.gz
#   /tmp/dv_hg002_wg_docker/intermediate/  (kept for forensics)

set -euo pipefail
cd "$(dirname "$0")/.."

DATA="${DV_GIAB_DIR:-/tmp/dv_giab}/full"
OUT="/tmp/dv_hg002_wg_docker"

if [ ! -f "${DATA}/HG002.bam" ] || [ ! -f "${DATA}/GRCh38.fa" ]; then
  echo "ERROR: missing data at ${DATA}" >&2
  exit 1
fi

if [ -f "${OUT}/output.vcf.gz" ]; then
  echo "==> ${OUT}/output.vcf.gz exists, skipping Docker run"
  exit 0
fi

mkdir -p "${OUT}"

echo "==> Docker DV 1.10.0 on HG002 WG (Rosetta 2, ~22 h ETA)"
date -u +"    Started: %Y-%m-%dT%H:%M:%SZ"

time docker run --rm \
  -v "${DATA}:/data:ro" \
  -v "${OUT}:/work" \
  google/deepvariant:1.10.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/GRCh38.fa \
    --reads=/data/HG002.bam \
    --output_vcf=/work/output.vcf.gz \
    --output_gvcf=/work/output.g.vcf.gz \
    --num_shards=14 \
    --intermediate_results_dir=/work/intermediate

date -u +"==> Done: %Y-%m-%dT%H:%M:%SZ"
echo "==> Docker WG VCF at ${OUT}/output.vcf.gz"
