#!/usr/bin/env bash
# Capture upstream DeepSomatic tumor-only reference outputs on the small
# chr20:10M-10.1M fixture for Docker-fidelity diffing.
#
# Modes covered:
#   WGS tumor-only      — deepsomatic.wgs_tumor_only model
#   FFPE_WGS tumor-only — deepsomatic.ffpe_wgs_tumor_only model
#
# Docker uses its internal PON at /opt/models/deepsomatic/pons/ automatically
# (via flags_for_calling in model.example_info.json). Our native binary uses
# the extracted PON at validation/work/deepsomatic_pon/ via --population_vcfs.
#
# Prerequisite: tools/reference/cache/ must contain:
#   HG002.chr20.10_10p1mb.bam{,.bai}
#   grch38_chr20.fasta{,.fai}
# (created by fetch_chr20_fixture.sh or capture_deeptrio_chr20.sh)
#
# Usage:
#   ./tools/reference/capture_deepsomatic_tumor_only_chr20.sh
#   DV_VERSION=1.10.0 REGION=chr20:10000000-10100000 \
#     ./tools/reference/capture_deepsomatic_tumor_only_chr20.sh
set -euo pipefail

cd "$(dirname "$0")"

DV_VERSION="${DV_VERSION:-1.10.0}"
REGION="${REGION:-chr20:10000000-10100000}"
IMG="google/deepsomatic:${DV_VERSION}"
OUT_BASE="output/deepsomatic_tumor_only"

# Resolve reference: prefer chr20-only fasta (faster), fall back to full GRCh38.
if [[ -s cache/grch38_chr20.fasta ]]; then
  REF_HOST="$(cd cache && pwd)"
  REF_BASE="grch38_chr20.fasta"
elif [[ -s /tmp/dv_giab/full/GRCh38.fa ]]; then
  REF_HOST="/tmp/dv_giab/full"
  REF_BASE="GRCh38.fa"
elif [[ -s /tmp/dv_giab/data/GRCh38.fa ]]; then
  REF_HOST="/tmp/dv_giab/data"
  REF_BASE="GRCh38.fa"
else
  echo "ERROR: no reference found — run fetch_chr20_fixture.sh first" >&2
  exit 1
fi

CACHE="$(cd cache && pwd)"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker not found" >&2
  exit 1
fi

echo "==> Pulling ${IMG} (linux/amd64)"
docker pull --platform linux/amd64 "${IMG}"

# ── Helper ───────────────────────────────────────────────────────────────────

run_deepsomatic_tumor_only() {
  local mode="$1"     # WGS_TUMOR_ONLY | FFPE_WGS_TUMOR_ONLY
  local mode_lc
  mode_lc="$(echo "${mode}" | tr '[:upper:]' '[:lower:]')"
  local out_dir="${OUT_BASE}/${mode_lc}"
  mkdir -p "${out_dir}"

  local out_abs
  out_abs="$(cd "${out_dir}" && pwd)"
  local done_flag="${out_abs}/docker.vcf.gz"
  if [[ -s "${done_flag}" ]]; then
    echo "==> SKIP (exists) ${mode} tumor-only Docker VCF"
    return 0
  fi

  echo "==> Running Docker deepsomatic ${mode} tumor-only on ${REGION}"

  docker run --rm \
    --platform linux/amd64 \
    -v "${CACHE}:/work/cache:ro" \
    -v "${REF_HOST}:/work/ref:ro" \
    -v "${out_abs}:/work/output" \
    "${IMG}" \
    bash -c "
      /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=${mode} \
        --ref=/work/ref/${REF_BASE} \
        --reads_tumor=/work/cache/HG002.chr20.10_10p1mb.bam \
        --output_vcf=/work/output/docker.vcf.gz \
        --regions=${REGION} \
        --intermediate_results_dir=/work/tmp_${mode_lc} \
        --num_shards=1
      tabix -f -p vcf /work/output/docker.vcf.gz
    "

  echo "==> Done: ${out_dir}/docker.vcf.gz"
}

# ── Run both modes ────────────────────────────────────────────────────────────

run_deepsomatic_tumor_only "WGS_TUMOR_ONLY"
run_deepsomatic_tumor_only "FFPE_WGS_TUMOR_ONLY"

echo ""
echo "==> Capture complete."
echo "  WGS tumor-only:      ${OUT_BASE}/wgs_tumor_only/docker.vcf.gz"
echo "  FFPE_WGS tumor-only: ${OUT_BASE}/ffpe_wgs_tumor_only/docker.vcf.gz"
