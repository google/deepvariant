#!/usr/bin/env bash
# Tier 1 — chr20 trio benchmark.
#
# Runs our native deepvariant on HG002, HG003, HG004 (chr20 only) and
# evaluates F1 vs GIAB v4.2.1 truth per sample via hap.py in Docker.
# All 3 samples use the same NovaSeq 35× PCR-free chr20 BAM family.
#
# Wall-time: ~30 min total on M4 Max (3 × 3 min run + 3 × 5 min hap.py).
#
# Outputs:
#   validation/output/<sample>_chr20/our.vcf.gz
#   validation/output/<sample>_chr20/happy.summary.csv
#   validation/output/chr20_trio_summary.tsv  (consolidated F1 table)

set -euo pipefail
cd "$(dirname "$0")/.."

DATA="${DV_GIAB_DIR:-/tmp/dv_giab/data}"
CKPT="${DV_CHECKPOINT:-/Users/benjamin/deepvariant/validation/work/wgs.dvw}"
SMALL="${DV_SMALL:-}"  # optional small-model path
INFER_BACKEND="${DV_INFERENCE:-metal}"

OUT_ROOT="validation/output"
SUMMARY="${OUT_ROOT}/chr20_trio_summary.tsv"
mkdir -p "${OUT_ROOT}"

run_sample() {
  local sample="$1"
  local truth_vcf="$2"
  local truth_bed="$3"
  local out="${OUT_ROOT}/${sample}_chr20"
  mkdir -p "${out}"

  echo
  echo "=========================================================="
  echo "==> ${sample} chr20 (NovaSeq 35× PCR-free)"
  echo "=========================================================="

  # Stage 1+2+3: deepvariant native run.
  if [ ! -f "${out}/our.vcf.gz" ]; then
    DV_ARGS=(
      --reads="${DATA}/${sample}.bam"
      --ref="${DATA}/GRCh38.fa"
      --regions=chr20
      --output_vcf="${out}/our.vcf.gz"
      --intermediate_results_dir="${out}/intermediate"
      --inference_backend="${INFER_BACKEND}"
      --model_type=WGS
      --checkpoint="${CKPT}"
      --num_shards=14
      --batch_size=512
    )
    [ -n "${SMALL}" ] && DV_ARGS+=(--small_model_path="${SMALL}")

    echo "==> ${sample}: native deepvariant on chr20"
    /usr/bin/time -h ./build-macos/bin/deepvariant run "${DV_ARGS[@]}" \
      2> "${out}/run_time.log" || {
        echo "ERROR: deepvariant run failed for ${sample}"; return 1; }
    grep -E "^real|^user|^sys" "${out}/run_time.log" || true
  else
    echo "==> ${out}/our.vcf.gz exists, skipping deepvariant run"
  fi

  # Index VCF for hap.py.
  if [ ! -f "${out}/our.vcf.gz.tbi" ]; then
    /opt/homebrew/bin/tabix -f -p vcf "${out}/our.vcf.gz"
  fi

  # Stage 2: hap.py vs GIAB v4.2.1 truth (chr20 subset via --location).
  # Note: explicit `--platform linux/amd64` triggers a Docker Desktop 500
  # error on this machine; omitting the flag lets Docker pick the image
  # platform (linux/amd64 by default) and emulate via Rosetta 2 / qemu.
  if [ ! -f "${out}/happy.summary.csv" ]; then
    echo "==> ${sample}: hap.py vs ${truth_vcf}"
    docker run --rm \
      -v "${DATA}:/data:ro" \
      -v "$(realpath "${out}"):/work" \
      jmcdani20/hap.py:v0.3.12 \
      /opt/hap.py/bin/hap.py \
        "/data/$(basename "${truth_vcf}")" \
        /work/our.vcf.gz \
        -f "/data/$(basename "${truth_bed}")" \
        -r /data/GRCh38.fa \
        -o /work/happy \
        --location chr20
  else
    echo "==> ${out}/happy.summary.csv exists, skipping hap.py"
  fi

  # Per-sample F1 summary.
  echo "==> ${sample}: F1 (PASS rows)"
  awk -F, 'NR==1 || /,PASS,/ {print}' "${out}/happy.summary.csv"
}

# HG002 — truth.vcf.gz / truth.bed (legacy filename in /tmp/dv_giab/data).
run_sample HG002 truth.vcf.gz       truth.bed
run_sample HG003 HG003.truth.vcf.gz HG003.truth.bed
run_sample HG004 HG004.truth.vcf.gz HG004.truth.bed

# Build consolidated F1 table.
echo
echo "=========================================================="
echo "==> Consolidated chr20 trio F1 table"
echo "=========================================================="
{
  printf "Sample\tType\tFilter\tTRUTH.TOTAL\tTRUTH.TP\tTRUTH.FN\tQUERY.TOTAL\tQUERY.FP\tRecall\tPrecision\tF1_Score\n"
  for s in HG002 HG003 HG004; do
    csv="${OUT_ROOT}/${s}_chr20/happy.summary.csv"
    if [ ! -f "${csv}" ]; then continue; fi
    awk -F, -v s="${s}" '
      NR > 1 && $2 == "PASS" {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               s, $1, $2, $3, $4, $5, $6, $7, $10, $11, $13
      }' "${csv}"
  done
} | tee "${SUMMARY}"

echo
echo "==> Summary saved at ${SUMMARY}"
