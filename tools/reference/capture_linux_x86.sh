#!/usr/bin/env bash
# Run upstream google/deepvariant Docker image under qemu linux/amd64 to
# capture reference outputs on a small chr20 region.
#
# Usage: ./capture_linux_x86.sh <wgs|wes|pacbio|ont_r104>
set -euo pipefail

cd "$(dirname "$0")"

VARIANT="${1:?usage: $0 <wgs|wes|pacbio|ont_r104>}"
DV_VERSION="1.10.0"
REGION="${REGION:-chr20:10000000-10100000}"
N_SHARDS="${N_SHARDS:-1}"

# Map variant -> upstream model_type flag
case "${VARIANT}" in
  wgs)       MODEL_TYPE="WGS" ;;
  wes)       MODEL_TYPE="WES" ;;
  pacbio)    MODEL_TYPE="PACBIO" ;;
  ont_r104)  MODEL_TYPE="ONT_R104" ;;
  *)         echo "unknown variant: ${VARIANT}" >&2; exit 2 ;;
esac

OUT="output/${VARIANT}"
mkdir -p "${OUT}"

if [[ ! -s cache/HG002.chr20.bam ]]; then
  echo "==> chr20 fixture not present; running fetch_chr20_fixture.sh"
  ./fetch_chr20_fixture.sh
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "error: docker not found" >&2
  exit 1
fi

# Make sure the linux/amd64 platform is available for qemu emulation.
docker buildx inspect default >/dev/null 2>&1 || docker buildx create --use --name dv-x86 || true

# Pull the upstream image (linux/amd64 explicit so qemu kicks in).
docker pull --platform linux/amd64 "google/deepvariant:${DV_VERSION}"

# We use the bundled run_deepvariant.py, but invoke each stage so we can capture
# intermediate TFRecords (examples + call_variants_outputs).
WORK="/work"
docker run --rm \
  --platform linux/amd64 \
  -v "${PWD}/cache:${WORK}/cache:ro" \
  -v "${PWD}/${OUT}:${WORK}/output" \
  "google/deepvariant:${DV_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type="${MODEL_TYPE}" \
    --ref="${WORK}/cache/grch38_chr20.fasta" \
    --reads="${WORK}/cache/HG002.chr20.bam" \
    --output_vcf="${WORK}/output/output.vcf.gz" \
    --output_gvcf="${WORK}/output/output.g.vcf.gz" \
    --num_shards="${N_SHARDS}" \
    --regions="${REGION}" \
    --intermediate_results_dir="${WORK}/output/intermediate" \
    2>&1 | tee "${OUT}/run.log"

# Extract the intermediate TFRecords for parity bench.
EXAMPLES_GLOB=("${OUT}"/intermediate/make_examples.tfrecord*)
CV_GLOB=("${OUT}"/intermediate/call_variants_output.tfrecord*)
if (( ${#EXAMPLES_GLOB[@]} > 0 && ${#CV_GLOB[@]} > 0 )); then
  cp "${EXAMPLES_GLOB[0]}" "${OUT}/examples_chr20.tfrecord"
  cp "${CV_GLOB[0]}"       "${OUT}/call_variants_chr20.tfrecord"
fi

[[ -s "${OUT}/examples_chr20.tfrecord" ]] || { echo "error: examples_chr20.tfrecord missing — make_examples produced no output" >&2; exit 1; }

# Build a 1000-example slice for fast bench iteration.
mkdir -p cache
python3 - <<PY "${OUT}/examples_chr20.tfrecord" "cache/${VARIANT}_chr20_1000.tfrecord"
import struct, sys
src, dst = sys.argv[1], sys.argv[2]
n_keep = 1000
with open(src, "rb") as fi, open(dst, "wb") as fo:
    n = 0
    while n < n_keep:
        ln_b = fi.read(8)
        if not ln_b: break
        (ln,) = struct.unpack("<Q", ln_b)
        crc1 = fi.read(4)
        payload = fi.read(ln)
        crc2 = fi.read(4)
        fo.write(ln_b); fo.write(crc1); fo.write(payload); fo.write(crc2)
        n += 1
print(f"wrote {n} examples to {dst}")
PY

# Manifest for traceability
{
  echo "{"
  echo "  \"variant\": \"${VARIANT}\","
  echo "  \"model_type\": \"${MODEL_TYPE}\","
  echo "  \"dv_version\": \"${DV_VERSION}\","
  echo "  \"region\": \"${REGION}\","
  echo "  \"captured_at\": \"$(date -u +%Y-%m-%dT%H:%M:%SZ)\","
  echo "  \"docker_platform\": \"linux/amd64 (qemu emulation)\""
  echo "}"
} > "${OUT}/manifest.json"

echo "==> done — ${OUT}"
ls -lh "${OUT}/"
