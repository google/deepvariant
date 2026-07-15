#!/usr/bin/env bash
# End-to-end HG002 WG comparison orchestrator.
#
# Stages:
#   1. Wait for our native WG to finish (validation/output/HG002_wg/our.vcf.gz)
#   2. Auto-launch Docker baseline (background, ~22 h)
#   3. When both available, compute:
#        - bit-level FILTER-class diff via diff_filter_classes.sh
#        - F1 ours vs GIAB v4.2.1 truth (already done by run_giab_wg_chunked.sh)
#        - F1 docker vs GIAB v4.2.1 truth (run hap.py separately on Docker output)
#        - PASS-set parity (bcftools isec)
#        - GT diff count on shared sites
#   4. Write consolidated report at validation/output/HG002_wg_benchmark.md
#
# Run this script BEFORE our native WG starts, or any time during. It
# polls until prerequisites are ready.

set -euo pipefail
cd "$(dirname "$0")/.."

OURS="validation/output/HG002_wg/our.vcf.gz"
DOCKER="/tmp/dv_hg002_wg_docker/output.vcf.gz"
TRUTH="/tmp/dv_giab/full/HG002.truth.vcf.gz"
TRUTH_BED="/tmp/dv_giab/full/HG002.truth.bed"
REF="/tmp/dv_giab/full/GRCh38.fa"
REPORT="validation/output/HG002_wg_benchmark.md"

mkdir -p "$(dirname "${REPORT}")"

echo "[$(date +%H:%M:%S)] Waiting for native WG (${OURS}) ..."
until [ -f "${OURS}" ] && [ -s "${OURS}" ]; do
  sleep 60
done
echo "[$(date +%H:%M:%S)] Native WG ready"

# Capture our native run's metadata.
NATIVE_PASS=$(zcat < "${OURS}" | awk '!/^#/ && $7=="PASS"' | wc -l)
NATIVE_TOTAL=$(zcat < "${OURS}" | grep -vc '^#')
NATIVE_F1_HAPPY="validation/output/HG002_wg/happy.summary.csv"

echo "[$(date +%H:%M:%S)] Native: ${NATIVE_TOTAL} sites, ${NATIVE_PASS} PASS"

# Stage 2: launch Docker baseline if not already done.
echo "[$(date +%H:%M:%S)] Launching Docker baseline (background, ~22 h) ..."
nohup ./validation/run_hg002_wg_docker_baseline.sh \
  > /tmp/hg002_wg_docker.log 2>&1 &
DOCKER_PID=$!
echo "[$(date +%H:%M:%S)] Docker baseline PID=${DOCKER_PID}"

echo "[$(date +%H:%M:%S)] Waiting for Docker WG (${DOCKER}) ..."
until [ -f "${DOCKER}" ] && [ -s "${DOCKER}" ]; do
  kill -0 "${DOCKER_PID}" 2>/dev/null || { echo "ERROR: Docker baseline (PID ${DOCKER_PID}) exited before producing ${DOCKER}; see /tmp/hg002_wg_docker.log" >&2; exit 1; }
  sleep 600   # 10 min poll — Docker takes ~22 h, no rush
done
echo "[$(date +%H:%M:%S)] Docker WG ready"

# Stage 3: comparisons.

# 3a. FILTER-class diff
DIFF_OUT=$(bash ./validation/diff_filter_classes.sh "${OURS}" "${DOCKER}" 2>&1)

# 3b. PASS-set parity via bcftools isec
ISEC_DIR="validation/output/HG002_wg/isec_vs_docker"
mkdir -p "${ISEC_DIR}"
/opt/homebrew/bin/bcftools isec -p "${ISEC_DIR}" \
  -Oz "${OURS}" "${DOCKER}"
SHARED=$(zcat < "${ISEC_DIR}/0002.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
ONLY_OURS=$(zcat < "${ISEC_DIR}/0000.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
ONLY_DOCKER=$(zcat < "${ISEC_DIR}/0001.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)

# 3c. PASS-set asymmetric diff
OURS_PASS_ONLY=$(zcat < "${ISEC_DIR}/0000.vcf.gz" 2>/dev/null | awk '!/^#/ && $7=="PASS"' | wc -l)
DOCKER_PASS_ONLY=$(zcat < "${ISEC_DIR}/0001.vcf.gz" 2>/dev/null | awk '!/^#/ && $7=="PASS"' | wc -l)

# 3d. GT diff on shared sites
GT_DIFF=$(paste \
  <(zcat < "${ISEC_DIR}/0002.vcf.gz" | grep -v '^#' | awk -F'\t' '{print $1"_"$2"_"$4"_"$5"\t"$10}' | cut -d':' -f1-2) \
  <(zcat < "${ISEC_DIR}/0003.vcf.gz" | grep -v '^#' | awk -F'\t' '{print $10}' | cut -d':' -f1) \
  | awk -F'\t' '{
      split($2, ours, ":"); split($3, doc, ":");
      if (ours[1] != doc[1]) c++
    } END {print c+0}')

# 3e. F1 hap.py for Docker output (parallel to our run's hap.py)
DOCKER_HAPPY="/tmp/dv_hg002_wg_docker/happy.summary.csv"
if [ ! -f "${DOCKER_HAPPY}" ]; then
  echo "[$(date +%H:%M:%S)] Running hap.py on Docker WG output ..."
  docker run --rm \
    -v /tmp/dv_giab/full:/data:ro \
    -v "$(realpath /tmp/dv_hg002_wg_docker):/work" \
    jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
      /data/HG002.truth.vcf.gz \
      /work/output.vcf.gz \
      -f /data/HG002.truth.bed \
      -r /data/GRCh38.fa \
      -o /work/happy
fi

# Stage 4: build report
cat > "${REPORT}" <<EOF
# HG002 whole-genome benchmark — Native arm64 vs Google Docker

**Date**: $(date -u +%Y-%m-%dT%H:%M:%SZ)
**Build commit**: $(git rev-parse --short HEAD)
**Hardware**: Apple M4 Max, 16 cores, 128 GB unified memory

## Site-set parity

| Metric | Count |
|--------|-------|
| Native total sites | ${NATIVE_TOTAL} |
| Docker total sites | $(zcat < "${DOCKER}" | grep -vc '^#') |
| Shared sites | ${SHARED} |
| Only ours | ${ONLY_OURS} |
| Only Docker | ${ONLY_DOCKER} |
| Native PASS | ${NATIVE_PASS} |
| Docker PASS | $(zcat < "${DOCKER}" | awk '!/^#/ && $7=="PASS"' | wc -l) |
| PASS-set ours-only | ${OURS_PASS_ONLY} |
| PASS-set docker-only | ${DOCKER_PASS_ONLY} |
| GT diffs on shared | ${GT_DIFF} |

## FILTER-class diff (on shared sites)

\`\`\`
${DIFF_OUT}
\`\`\`

## F1 vs GIAB v4.2.1 truth

### Native

\`\`\`
$(awk -F, 'NR==1 || /,PASS,/' "${NATIVE_F1_HAPPY}" 2>&1)
\`\`\`

### Docker

\`\`\`
$(awk -F, 'NR==1 || /,PASS,/' "${DOCKER_HAPPY}" 2>&1)
\`\`\`

## Verdict

Generated $(date -u +%Y-%m-%dT%H:%M:%SZ).
EOF

echo
echo "==> Report at ${REPORT}"
cat "${REPORT}"
