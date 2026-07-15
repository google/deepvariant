#!/usr/bin/env bash
# Pull the small GIAB chr20 fixture used by the reference capture.
# This is the same dataset Google uses in the public DeepVariant case study.
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p cache

REF_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/grch38_chr20.fasta"
REF_FAI_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/grch38_chr20.fasta.fai"
BAM_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.chr20.bam"
BAI_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.chr20.bam.bai"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "${out}" ]]; then
    echo "  ${out} present, skip"
    return
  fi
  echo "  fetching ${url}"
  curl -fL --retry 3 --connect-timeout 15 -o "${out}.partial" "${url}"
  mv "${out}.partial" "${out}"
}

echo "==> chr20 reference fixture"
fetch "${REF_URL}"     cache/grch38_chr20.fasta
fetch "${REF_FAI_URL}" cache/grch38_chr20.fasta.fai
fetch "${BAM_URL}"     cache/HG002.chr20.bam
fetch "${BAI_URL}"     cache/HG002.chr20.bam.bai

echo "==> done"
ls -lh cache/
