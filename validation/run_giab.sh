#!/usr/bin/env bash
# Run our deepvariant binary on a GIAB benchmark sample and compute F1
# scores via hap.py. Targets the plan's scientific gates:
#   SNP F1   ≥ upstream F1 − 0.05 %
#   INDEL F1 ≥ upstream F1 − 0.10 %
#
# hap.py runs in Docker (it's a validation tool, not part of the binary),
# so this script needs Docker. It does NOT run our deepvariant in Docker.
#
# Usage:
#   ./validation/run_giab.sh <sample-id> [region|""]
#   ./validation/run_giab.sh HG002 chr20            # chr20 only (~3-4 min)
#   ./validation/run_giab.sh HG002 ""               # whole-genome (~3 hours)
#   DV_VALIDATION_OUT=/tmp/wg_hg002 ./validation/run_giab.sh HG002 ""
#
# Phase 9 / Step 5 — empty `region` arg invokes whole-genome mode (no
# --regions filter passed to deepvariant + no --location filter to
# hap.py). Wall-time on M4 Max ~14-thread: ~50× chr20 = ~3 hours per
# sample. F1 evaluated over the full GIAB high-confidence regions
# (truth.bed is already whole-genome).
#
# Inputs expected at $DV_GIAB_DIR (default /tmp/dv_giab/data):
#   <sample>.bam{,.bai}                — aligned reads
#   GRCh38.fa{,.fai}                   — reference (full or chr20-only)
#   truth.vcf.gz{,.tbi}                — GIAB benchmark VCF
#   truth.bed                          — GIAB high-confidence regions

set -euo pipefail
cd "$(dirname "$0")/.."

SAMPLE="${1:?usage: $0 <sample> [region|\"\" for whole-genome]}"
# Default region is chr20 unless explicitly passed (empty string allowed
# for whole-genome).
if [ "$#" -ge 2 ]; then
  REGION="$2"
else
  REGION="chr20"
fi
DATA="${DV_GIAB_DIR:-/tmp/dv_giab/data}"
OUT="${DV_VALIDATION_OUT:-validation/output/${SAMPLE}}"

mkdir -p "${OUT}"
if [ -z "${REGION}" ]; then
  echo "==> Running our deepvariant on ${SAMPLE} WHOLE-GENOME (~3 h on M4 Max)"
else
  echo "==> Running our deepvariant on ${SAMPLE} ${REGION}"
fi

# Build deepvariant args: omit --regions when whole-genome mode.
DV_ARGS=(
  --reads="${DATA}/${SAMPLE}.bam"
  --ref="${DATA}/GRCh38.fa"
  --output_vcf="${OUT}/our.vcf"
  --num_shards=1
  --intermediate_results_dir="${OUT}/intermediate"
  --model=tools/conversion/models/wgs.mlpackage
  --small_model_path=tools/conversion/models/wgs_small.mlpackage
  --compute_units=all
)
if [ -n "${REGION}" ]; then
  DV_ARGS+=(--regions="${REGION}")
fi
time ./build-macos/bin/deepvariant run "${DV_ARGS[@]}"

# Compress + index for hap.py
"${HOMEBREW_PREFIX:-/opt/homebrew}/bin/bgzip" -f "${OUT}/our.vcf"
"${HOMEBREW_PREFIX:-/opt/homebrew}/bin/tabix" -f -p vcf "${OUT}/our.vcf.gz"

echo "==> hap.py vs GIAB truth (Docker linux/amd64)"
HAPPY_ARGS=(
  /data/truth.vcf.gz
  /work/our.vcf.gz
  -f /data/truth.bed
  -r /data/GRCh38.fa
  -o /work/happy
)
if [ -n "${REGION}" ]; then
  HAPPY_ARGS+=(--location "${REGION}")
fi
docker run --rm --platform linux/amd64 \
  -v "${DATA}:/data:ro" \
  -v "$(realpath "${OUT}"):/work" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py "${HAPPY_ARGS[@]}"

echo
echo "==> F1 scores"
"${HOMEBREW_PREFIX:-/opt/homebrew}/bin/csvcut" -c "Type,Filter,METRIC.Recall,METRIC.Precision,METRIC.F1_Score" \
  "${OUT}/happy.summary.csv" 2>/dev/null \
  || cat "${OUT}/happy.summary.csv"
