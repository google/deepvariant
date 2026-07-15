#!/usr/bin/env bash
# Phase 8 / Tier 2 — Multi-seed Test-Time Augmentation orchestrator.
#
# Runs `deepvariant run` N times with different --tta_seed_offset values,
# each producing a slightly different pileup view of the same reads
# (different shuffle in DownsampleReadIndices when coverage > pileup
# height, different reservoir sample). Then aggregates the N output
# VCFs by majority vote on FILTER class per site.
#
# Expected effect: +0.05-0.20 % F1 lift on borderline GQ ≈ 20 / QUAL ≈ 1
# sites (Shorten & Khoshgoftaar 2019 J. Big Data — typical TTA gain in
# CV literature; never published for variant calling specifically but
# the mechanism applies). Cost: N× wall-time (~N×4 min/chr20).
#
# Trade-off: TTA averaging may introduce its own noise at boundaries.
# If majority vote disagrees with the baseline call, the TTA call wins.
# Net F1 gain depends on whether borderline sites are true-positive or
# true-negative — empirical measurement required.
#
# Usage:
#   ./validation/run_tta.sh <bam> <ref> <out_dir> <region> [N=5]
#
# Outputs in <out_dir>/:
#   pass-K/ (K = 0..N-1)             — per-pass deepvariant output
#   tta_majority.vcf.gz              — majority-vote final VCF
#   tta_summary.tsv                  — per-site disagreement count

set -euo pipefail
cd "$(dirname "$0")/.."

if [ "$#" -lt 4 ]; then
  echo "usage: $0 <bam> <ref> <out_dir> <region> [N=5]" >&2
  exit 2
fi

BAM="$1"
REF="$2"
OUT="$3"
REGION="$4"
N="${5:-5}"

DV="${DV_BIN:-./build-macos/bin/deepvariant}"
MODEL="${DV_MODEL:-validation/work/wgs.dvw}"
BCFTOOLS="${BCFTOOLS:-/opt/homebrew/bin/bcftools}"

mkdir -p "${OUT}"

echo "==> Multi-seed TTA on ${REGION} with N=${N} passes"
echo "    BAM:    ${BAM}"
echo "    REF:    ${REF}"
echo "    MODEL:  ${MODEL}"
echo

for ((K = 0; K < N; ++K)); do
  PASS_DIR="${OUT}/pass-${K}"
  if [ -f "${PASS_DIR}/output.vcf.gz" ]; then
    echo "==> Pass ${K}: cached at ${PASS_DIR}, skipping"
    continue
  fi
  mkdir -p "${PASS_DIR}"
  echo "==> Pass ${K} (tta_seed_offset=${K})"
  time "${DV}" run \
    --reads="${BAM}" \
    --ref="${REF}" \
    --output_vcf="${PASS_DIR}/output.vcf.gz" \
    --model_type=WGS \
    --model="${MODEL}" \
    --regions="${REGION}" \
    --tta_seed_offset="${K}" \
    --intermediate_results_dir="${PASS_DIR}/intermediate" \
    > "${PASS_DIR}/run.log" 2>&1
  "${BCFTOOLS}" index -t -f "${PASS_DIR}/output.vcf.gz"
done

echo
echo "==> Aggregating ${N} passes by majority-vote per site"

# Build a per-pass FILTER table: CHROM\tPOS\tFILTER per pass.
for ((K = 0; K < N; ++K)); do
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%FILTER\n' \
    "${OUT}/pass-${K}/output.vcf.gz" > "${OUT}/pass-${K}.tsv"
done

# Pick majority FILTER per (CHROM,POS); ties go to baseline (pass 0).
python3 - <<PYEOF
import collections, csv, gzip, os, sys

n = ${N}
out_dir = "${OUT}"

# Collect per-site FILTERs across passes.
agg = collections.defaultdict(list)
for k in range(n):
    with open(f"{out_dir}/pass-{k}.tsv") as f:
        for line in f:
            chrom, pos, flt = line.rstrip("\n").split("\t")
            agg[(chrom, int(pos))].append(flt)

# Majority vote per site.
disagree = 0
votes_sum = collections.Counter()
with open(f"{out_dir}/tta_summary.tsv", "w") as fout:
    fout.write("CHROM\tPOS\tBASELINE\tWINNER\tVOTES\n")
    for key in sorted(agg.keys()):
        flts = agg[key]
        c = collections.Counter(flts)
        baseline = flts[0]
        winner, _ = c.most_common(1)[0]
        votes_sum[winner] += 1
        if winner != baseline:
            disagree += 1
        fout.write(f"{key[0]}\t{key[1]}\t{baseline}\t{winner}\t"
                   f"{','.join(f'{f}={v}' for f, v in c.most_common())}\n")

print(f"  total sites:     {len(agg)}")
print(f"  disagreement:    {disagree} ({100*disagree/max(1,len(agg)):.3f} %)")
print(f"  winners by FILTER: {dict(votes_sum)}")
PYEOF

echo
echo "==> Per-site summary: ${OUT}/tta_summary.tsv"
echo "==> Pass-0 baseline: ${OUT}/pass-0/output.vcf.gz"
echo
echo "Note: this script reports majority-vote stats only. To produce a"
echo "merged VCF where FILTER is replaced by the majority-vote winner,"
echo "extend this script with a bcftools annotate pass — TODO once the"
echo "F1 effect of majority voting has been measured."
