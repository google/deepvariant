#!/usr/bin/env bash
# Diff two VCFs by per-site FILTER class (PASS / RefCall / NoCall / LowQual).
#
# Standard tool used throughout this project to verify Docker FILTER-class
# parity. Outputs:
#   - shared site count
#   - site-set divergence (only_a / only_b)
#   - FM total on shared sites
#   - per-transition histogram (e.g., "PASS -> NoCall: 32")
#
# Usage:
#   ./validation/diff_filter_classes.sh <ours.vcf.gz> <docker.vcf.gz>
#
# Optional environment:
#   BCFTOOLS=/path/to/bcftools  (default: /opt/homebrew/bin/bcftools)

set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "usage: $0 <ours.vcf.gz> <docker.vcf.gz>" >&2
  exit 1
fi

OURS="$1"
DOCKER="$2"
BCFTOOLS="${BCFTOOLS:-/opt/homebrew/bin/bcftools}"
WORK="$(mktemp -d)"
trap '[[ -n "${DV_DIFF_KEEP:-}" ]] || rm -rf "${WORK}"' EXIT

# Index inputs if not already.
[ -f "${OURS}.tbi" ]   || "${BCFTOOLS}" index -t -f "${OURS}"
[ -f "${DOCKER}.tbi" ] || "${BCFTOOLS}" index -t -f "${DOCKER}"

"${BCFTOOLS}" isec -p "${WORK}" -Ov "${OURS}" "${DOCKER}" 2>&1 | tail -2

ONLY_OURS=$(grep -cv "^#" "${WORK}/0000.vcf" || true)
ONLY_DOCKER=$(grep -cv "^#" "${WORK}/0001.vcf" || true)
SHARED=$(grep -cv "^#" "${WORK}/0002.vcf" || true)

"${BCFTOOLS}" query -f '%CHROM\t%POS\t%FILTER\n' "${WORK}/0002.vcf" > "${WORK}/ours_filt.tsv"
"${BCFTOOLS}" query -f '%CHROM\t%POS\t%FILTER\n' "${WORK}/0003.vcf" > "${WORK}/docker_filt.tsv"

FM=$(paste "${WORK}/ours_filt.tsv" "${WORK}/docker_filt.tsv" | awk -F'\t' '$3 != $6' | wc -l | tr -d ' ')

echo
echo "=== FILTER-class diff: ${OURS} vs ${DOCKER} ==="
printf '  shared sites    : %s\n' "${SHARED}"
printf '  only ours       : %s\n' "${ONLY_OURS}"
printf '  only docker     : %s\n' "${ONLY_DOCKER}"
printf '  FM on shared    : %s\n' "${FM}"

if [ "${FM}" -gt 0 ]; then
  echo
  echo "=== transition histogram ==="
  paste "${WORK}/ours_filt.tsv" "${WORK}/docker_filt.tsv" |
    awk -F'\t' '$3 != $6 {print "  "$3" -> "$6}' |
    sort | uniq -c | sort -rn
fi

# 0 FM + 0 site-set diff = bit-identical FILTER classes (the gate).
if [ "${ONLY_OURS}" -eq 0 ] && [ "${ONLY_DOCKER}" -eq 0 ] && [ "${FM}" -eq 0 ]; then
  echo
  echo "✅ 100 % FILTER-class parity"
  exit 0
fi

if [ -n "${DV_DIFF_KEEP:-}" ]; then
  echo
  echo "Intermediate files kept at: ${WORK}"
fi

exit 1
