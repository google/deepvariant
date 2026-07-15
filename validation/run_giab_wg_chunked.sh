#!/usr/bin/env bash
# Tier 2 — whole-genome trio benchmark, chunked execution.
#
# Disk-managed pipeline that processes one chromosome at a time per
# sample, freeing intermediate files (examples.tfrecord ~50 GB) between
# chunks. Required because chr20 alone makes 12 GB examples; full WG
# would make ~600 GB without chunking, blowing our 127 GB disk budget.
#
# Wall-time per sample: ~2.5 h compute + ~30 min hap.py = ~3 h
# Sequential 3 samples: ~9 h.
#
# Disk peak per sample:
#   BAM (40 GB) + chr1 examples (~50 GB) + accumulated VCF (~200 MB)
#   ≈ 90 GB peak (fits in 127 GB after subtract reference 3 GB).
#
# Usage:
#   ./validation/run_giab_wg_chunked.sh                  # all 3 samples
#   ./validation/run_giab_wg_chunked.sh HG002            # single sample
#   DV_GIAB_DIR=/tmp/dv_giab/full DV_KEEP_BAM=1 ./validation/run_giab_wg_chunked.sh
#
# Prereq: ./validation/download_giab_full_genome.sh ran (~120 GB at
#   ${DV_GIAB_DIR}/full).

set -euo pipefail
cd "$(dirname "$0")/.."

DATA="${DV_GIAB_DIR:-/tmp/dv_giab}/full"
CKPT="${DV_CHECKPOINT:-/Users/benjamin/deepvariant/validation/work/wgs.dvw}"
SMALL="${DV_SMALL:-/Users/benjamin/deepvariant/validation/work/wgs_small_weights}"
INFER_BACKEND="${DV_INFERENCE:-metal}"
# When INFER_BACKEND=ane_speculate, --checkpoint must be the .mlpackage and
# --ane_speculate_metal_checkpoint must be the .dvw rerun bundle.
ANE_DVW="${DV_ANE_DVW:-/Users/benjamin/deepvariant/validation/work/wgs.dvw}"
ANE_CONF="${DV_ANE_CONF:-0.99}"
KEEP_BAM="${DV_KEEP_BAM:-0}"          # set 1 to retain BAM after sample done
NUM_SHARDS="${DV_NUM_SHARDS:-14}"   # M4 Max has 14 P-cores; saturate make_examples
BATCH_SIZE="${DV_BATCH_SIZE:-512}"

if [ ! -f "${DATA}/GRCh38.fa" ]; then
  echo "ERROR: missing reference at ${DATA}/GRCh38.fa" >&2
  echo "       Run ./validation/download_giab_full_genome.sh first." >&2
  exit 1
fi

# Chunks: 25 chunks (chr1-22, chrX, chrY, chrM). Order matters only for
# the intermediate disk peak — chr1 is the largest, do it first when
# the BAM was just downloaded (no other intermediates lying around).
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
        chrX chrY chrM)

run_chunk() {
  local sample="$1" chrom="$2" out_dir="$3"
  local chunk_vcf="${out_dir}/chunks/${chrom}.vcf.gz"
  local inter_dir="${out_dir}/intermediate_${chrom}"

  if [ -f "${chunk_vcf}" ] && [ -s "${chunk_vcf}" ]; then
    echo "    ${sample}/${chrom}: cached, skipping"
    return 0
  fi

  mkdir -p "$(dirname "${chunk_vcf}")" "${inter_dir}"

  echo "    ${sample}/${chrom}: deepvariant run …"
  local sm_args=()
  [ -n "${SMALL}" ] && [ -d "${SMALL}" ] && sm_args+=(--small_model_path="${SMALL}")
  local ane_args=()
  if [ "${INFER_BACKEND}" = "ane_speculate" ]; then
    ane_args+=(--ane_speculate_metal_checkpoint="${ANE_DVW}")
    ane_args+=(--ane_speculate_confidence="${ANE_CONF}")
  fi
  /usr/bin/time -p ./build-macos/bin/deepvariant run \
    --reads="${DATA}/${sample}.bam" \
    --ref="${DATA}/GRCh38.fa" \
    --regions="${chrom}" \
    --output_vcf="${chunk_vcf}" \
    --intermediate_results_dir="${inter_dir}" \
    --inference_backend="${INFER_BACKEND}" \
    --model_type=WGS \
    --checkpoint="${CKPT}" \
    --num_shards="${NUM_SHARDS}" \
    --batch_size="${BATCH_SIZE}" \
    "${sm_args[@]+"${sm_args[@]}"}" \
    "${ane_args[@]+"${ane_args[@]}"}" \
    > "${inter_dir}/run.log" 2>&1

  # Free disk: drop intermediate examples.tfrecord (the 50 GB beast)
  # but keep run.log + cvo.tfrecord for forensics if requested.
  if [ "${DV_KEEP_INTERMEDIATE:-0}" != "1" ]; then
    rm -rf "${inter_dir}"
  fi

  /opt/homebrew/bin/tabix -f -p vcf "${chunk_vcf}"
}

run_sample() {
  local sample="$1" truth_vcf="$2" truth_bed="$3"
  local out="validation/output/${sample}_wg${DV_OUT_SUFFIX:-}"
  mkdir -p "${out}/chunks"

  echo
  echo "============================================================"
  echo "==> ${sample} whole-genome (chunked)"
  echo "============================================================"
  echo "    BAM:   ${DATA}/${sample}.bam"
  echo "    Truth: ${truth_vcf}"
  date -u +"    Started: %Y-%m-%dT%H:%M:%SZ"
  echo

  if [ ! -f "${DATA}/${sample}.bam" ]; then
    echo "ERROR: missing ${DATA}/${sample}.bam" >&2
    return 1
  fi

  local sample_start sample_end
  sample_start="$(date +%s)"

  for chrom in "${CHROMS[@]}"; do
    run_chunk "${sample}" "${chrom}" "${out}"
    df -g "${out}" 2>/dev/null \
      | awk -v c="${chrom}" 'NR==2{printf "    [disk after %s] %d GB free\n", c, $4}'
  done

  # Concatenate chunk VCFs in chrom order.
  if [ ! -f "${out}/our.vcf.gz" ]; then
    echo "    Concatenating ${#CHROMS[@]} chunk VCFs …"
    /opt/homebrew/bin/bcftools concat -Oz \
      $(for c in "${CHROMS[@]}"; do
          [ -f "${out}/chunks/${c}.vcf.gz" ] && echo "${out}/chunks/${c}.vcf.gz"
        done) \
      > "${out}/our.vcf.gz"
    /opt/homebrew/bin/tabix -f -p vcf "${out}/our.vcf.gz"
  fi

  sample_end="$(date +%s)"
  local elapsed=$((sample_end - sample_start))
  echo "    deepvariant wall-time: ${elapsed} s ($((elapsed / 60)) min)"

  # hap.py vs GIAB v4.2.1 truth, whole-genome (no --location).
  # Note: omit explicit `--platform linux/amd64` — that flag triggers a
  # Docker Desktop 500 error; default platform selection works.
  if [ ! -f "${out}/happy.summary.csv" ]; then
    echo "    hap.py vs ${truth_vcf} (whole-genome) …"
    docker run --rm \
      -v "${DATA}:/data:ro" \
      -v "$(realpath "${out}"):/work" \
      jmcdani20/hap.py:v0.3.12 \
      /opt/hap.py/bin/hap.py \
        "/data/$(basename "${truth_vcf}")" \
        /work/our.vcf.gz \
        -f "/data/$(basename "${truth_bed}")" \
        -r /data/GRCh38.fa \
        -o /work/happy
  fi

  echo
  echo "    F1 (PASS rows):"
  awk -F, 'NR==1 || /,PASS,/ {print "      "$0}' "${out}/happy.summary.csv"

  # Optional: free up the BAM after the sample completes.
  if [ "${KEEP_BAM}" != "1" ]; then
    echo "    Removing ${DATA}/${sample}.bam (set DV_KEEP_BAM=1 to retain)"
    rm -f "${DATA}/${sample}.bam" "${DATA}/${sample}.bam.bai"
  fi
}

# Single-sample mode (passed as first arg) or all 3.
if [ "$#" -ge 1 ]; then
  case "$1" in
    HG002) run_sample HG002 "${DATA}/HG002.truth.vcf.gz" "${DATA}/HG002.truth.bed" ;;
    HG003) run_sample HG003 "${DATA}/HG003.truth.vcf.gz" "${DATA}/HG003.truth.bed" ;;
    HG004) run_sample HG004 "${DATA}/HG004.truth.vcf.gz" "${DATA}/HG004.truth.bed" ;;
    *) echo "Unknown sample: $1 (expect HG002 / HG003 / HG004)"; exit 1 ;;
  esac
else
  run_sample HG002 "${DATA}/HG002.truth.vcf.gz" "${DATA}/HG002.truth.bed"
  run_sample HG003 "${DATA}/HG003.truth.vcf.gz" "${DATA}/HG003.truth.bed"
  run_sample HG004 "${DATA}/HG004.truth.vcf.gz" "${DATA}/HG004.truth.bed"
fi

echo
echo "============================================================"
echo "==> Trio whole-genome F1 validation complete."
echo "    Results at validation/output/{HG002,HG003,HG004}_wg/"
echo "============================================================"

# Build consolidated F1 table.
SUMMARY="validation/output/wg_trio_summary.tsv"
{
  printf "Sample\tType\tFilter\tTRUTH.TOTAL\tTRUTH.TP\tTRUTH.FN\tQUERY.TOTAL\tQUERY.FP\tRecall\tPrecision\tF1_Score\n"
  for s in HG002 HG003 HG004; do
    csv="validation/output/${s}_wg/happy.summary.csv"
    [ ! -f "${csv}" ] && continue
    awk -F, -v s="${s}" '
      NR > 1 && $2 == "PASS" {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               s, $1, $2, $3, $4, $5, $6, $7, $10, $11, $13
      }' "${csv}"
  done
} | tee "${SUMMARY}"

echo
echo "==> Consolidated WG F1 saved at ${SUMMARY}"
