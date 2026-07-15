#!/usr/bin/env bash
# Capture upstream DeepTrio reference outputs on a small chr20 fixture for
# Docker-fidelity diffing during the Step 1 (DeepTrio) port. Mirrors the
# trio quick-start at docs/deeptrio-quick-start.md but kept minimal so it
# runs in ~10 min under qemu linux/amd64 on Apple Silicon.
#
# Outputs into tools/reference/cache/ (BAMs/FASTA, shared with the WGS
# capture) and tools/reference/output/deeptrio/ (per-sample VCFs +
# intermediate TFRecords for stage-by-stage diffing).
#
# Usage: ./capture_deeptrio_chr20.sh
set -euo pipefail

cd "$(dirname "$0")"

DV_VERSION="${DV_VERSION:-1.10.0}"
REGION="${REGION:-chr20:10000000-10100000}"
N_SHARDS="${N_SHARDS:-1}"
OUT="output/deeptrio"
mkdir -p "${OUT}" cache

HTTPDIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"
REF_FTP="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"

fetch() {
  local url="$1" out="$2"
  if [[ -s "${out}" ]]; then
    echo "  ${out} present, skip"
    return
  fi
  echo "  fetching ${url}"
  /usr/bin/curl -fL --retry 3 --connect-timeout 15 -o "${out}.partial" "${url}"
  mv "${out}.partial" "${out}"
}

echo "==> chr20 trio fixture"
for s in HG002 HG003 HG004; do
  fetch "${HTTPDIR}/${s}.chr20.10_10p1mb.bam"     "cache/${s}.chr20.10_10p1mb.bam"
  fetch "${HTTPDIR}/${s}.chr20.10_10p1mb.bam.bai" "cache/${s}.chr20.10_10p1mb.bam.bai"
done
# Reuse the WGS fixture's reference (deeptrio uses full GRCh38 but for
# our tiny region either works; if grch38_chr20.fasta is present, use it,
# otherwise fall back to the full GRCh38 already on disk for WGS).
if [[ -s cache/grch38_chr20.fasta ]]; then
  REF=cache/grch38_chr20.fasta
elif [[ -s /tmp/dv_giab/data/GRCh38.fa ]]; then
  REF=/tmp/dv_giab/data/GRCh38.fa
else
  echo "  fetching full GRCh38 (slow, but only once)"
  fetch "${REF_FTP}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" cache/GRCh38.fa.fai
  echo "  GRCh38.fasta itself not auto-fetched here — set DV_REF env var."
  exit 1
fi
echo "  reference: ${REF}"

if ! command -v docker >/dev/null 2>&1; then
  echo "error: docker not found" >&2
  exit 1
fi

# linux/amd64 needed; pull (qemu emulation under Apple Silicon)
docker pull --platform linux/amd64 "google/deepvariant:deeptrio-${DV_VERSION}"

# Set up the docker mount roots; absolute path to the reference too.
WORK="/work"
REF_HOST="$(cd "$(dirname "${REF}")" && pwd)"
REF_BASE="$(basename "${REF}")"

docker run --rm \
  --platform linux/amd64 \
  -v "${PWD}/cache:${WORK}/cache:ro" \
  -v "${REF_HOST}:${WORK}/ref:ro" \
  -v "${PWD}/${OUT}:${WORK}/output" \
  "google/deepvariant:deeptrio-${DV_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
    --model_type=WGS \
    --ref="${WORK}/ref/${REF_BASE}" \
    --reads_child="${WORK}/cache/HG002.chr20.10_10p1mb.bam" \
    --reads_parent1="${WORK}/cache/HG003.chr20.10_10p1mb.bam" \
    --reads_parent2="${WORK}/cache/HG004.chr20.10_10p1mb.bam" \
    --output_vcf_child="${WORK}/output/HG002.output.vcf.gz" \
    --output_vcf_parent1="${WORK}/output/HG003.output.vcf.gz" \
    --output_vcf_parent2="${WORK}/output/HG004.output.vcf.gz" \
    --output_gvcf_child="${WORK}/output/HG002.g.vcf.gz" \
    --output_gvcf_parent1="${WORK}/output/HG003.g.vcf.gz" \
    --output_gvcf_parent2="${WORK}/output/HG004.g.vcf.gz" \
    --sample_name_child=HG002 \
    --sample_name_parent1=HG003 \
    --sample_name_parent2=HG004 \
    --num_shards="${N_SHARDS}" \
    --regions="${REGION}" \
    --intermediate_results_dir="${WORK}/output/intermediate" \
    2>&1 | tee "${OUT}/run.log"

echo "==> done. Outputs in ${OUT}/"
ls -la "${OUT}/" | head -20
