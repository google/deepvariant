#!/usr/bin/env bash
# Extract .dvw weight bundles and small-model .npy files for every
# DeepVariant / DeepTrio / DeepSomatic model variant.
#
# The .dvw format is the runtime weight file consumed by the Phase 5.5
# Metal/BNNS inference path (see tools/conversion/extract_weights.py).
# Small-model .npy files are consumed by the BNNS-CPU FP32 MLP
# (see tools/conversion/extract_small_model_weights.sh).
#
# ── Big-model (.dvw) strategy ─────────────────────────────────────────────
# SavedModels are already fetched into tools/conversion/models/<name>/ by
# fetch_savedmodel.sh (evidenced by all .mlpackage conversions already done).
# extract_weights.py is TF-free (pure-protobuf via tensor_bundle_reader.py),
# so it runs under any Python 3 that has numpy — we reuse the DeepVariant
# Docker (google/deepvariant:1.10.0) which already has numpy 1.24.
#
# ── Small-model (.npy) strategy ───────────────────────────────────────────
# Small models live ONLY inside Docker images; they are not on GCS.
# We reuse tools/conversion/extract_small_model_weights.sh which spins up
# the appropriate Docker image and exports Keras layer weights as .npy files.
#
# Docker images:
#   google/deepvariant:1.10.0          — WGS, WES, PacBio, ONT, Hybrid,
#                                        MaSeq, RNASeq big+small models
#   google/deepvariant:deeptrio-1.10.0 — DeepTrio big+small models
#   google/deepsomatic:1.10.0          — DeepSomatic big+small models
#
# ── Outputs ───────────────────────────────────────────────────────────────
# validation/work/<name>.dvw              big-model weight bundles
# validation/work/<name>_small/           small-model .npy directories
#
# ── Idempotency ───────────────────────────────────────────────────────────
# Each target is skipped if the output file/directory already exists and is
# non-empty. Re-run freely; partial runs resume from the first missing file.
#
# Usage:
#   tools/reference/extract_all_model_weights.sh
#
# Override output directory:
#   DVW_OUT=/path/to/dir tools/reference/extract_all_model_weights.sh
#
# Override DV version (default 1.10.0):
#   DV_VERSION=1.10.0 tools/reference/extract_all_model_weights.sh

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MODELS_DIR="${REPO_ROOT}/tools/conversion/models"
CONVERSION_DIR="${REPO_ROOT}/tools/conversion"
DVW_OUT="${DVW_OUT:-${REPO_ROOT}/validation/work}"
DV_VERSION="${DV_VERSION:-1.10.0}"

EXTRACT_SMALL="${CONVERSION_DIR}/extract_small_model_weights.sh"

mkdir -p "${DVW_OUT}"

# ── Docker images ─────────────────────────────────────────────────────────
IMG_DV="google/deepvariant:${DV_VERSION}"
IMG_DT="google/deepvariant:deeptrio-${DV_VERSION}"
IMG_DS="google/deepsomatic:${DV_VERSION}"

# ── Helpers ───────────────────────────────────────────────────────────────

log() { echo "==> $*"; }

# extract_dvw <savedmodel_dir> <out.dvw> [docker_image]
# Runs extract_weights.py inside Docker so we get numpy 1.24 + the
# vendored tensor_bundle_reader without polluting the host venv.
# The script is TF-free; any Docker image with numpy works — we use IMG_DV.
extract_dvw() {
  local src_dir="$1"   # absolute host path to SavedModel directory
  local out_dvw="$2"   # absolute host path for output .dvw
  local docker_img="${3:-${IMG_DV}}"

  if [[ -s "${out_dvw}" ]]; then
    log "SKIP (exists) ${out_dvw##*/}"
    return 0
  fi

  log "Extracting $(basename "${out_dvw}") from $(basename "${src_dir}")"

  local out_dir
  out_dir="$(dirname "${out_dvw}")"
  local out_name
  out_name="$(basename "${out_dvw}")"

  docker run --rm --platform linux/amd64 \
    -v "${CONVERSION_DIR}:/conversion:ro" \
    -v "${src_dir}:/savedmodel:ro" \
    -v "${out_dir}:/out" \
    "${docker_img}" \
    python3 /conversion/extract_weights.py /savedmodel "/out/${out_name}"
}

# extract_small <out_dir> <model_path_in_docker> <docker_tag>
# Delegates to extract_small_model_weights.sh (uses google/deepvariant:<tag>).
extract_small() {
  local out_dir="$1"         # host path for .npy output directory
  local model_path="$2"      # path inside Docker, e.g. /opt/smallmodels/wgs/model.keras
  local docker_tag="$3"      # third arg to extract_small_model_weights.sh

  if [[ -d "${out_dir}" && -n "$(ls -A "${out_dir}" 2>/dev/null)" ]]; then
    log "SKIP (exists) $(basename "${out_dir}")"
    return 0
  fi

  log "Extracting small model → $(basename "${out_dir}") (${model_path})"
  bash "${EXTRACT_SMALL}" "${out_dir}" "${model_path}" "${docker_tag}"
}

# extract_small_from_image <out_dir> <model_path_in_docker> <docker_image>
# Like extract_small but uses an arbitrary Docker image directly (for
# DeepSomatic which uses google/deepsomatic, not google/deepvariant).
extract_small_from_image() {
  local out_dir="$1"
  local model_path="$2"
  local docker_img="$3"

  if [[ -d "${out_dir}" && -n "$(ls -A "${out_dir}" 2>/dev/null)" ]]; then
    log "SKIP (exists) $(basename "${out_dir}")"
    return 0
  fi

  log "Extracting small model → $(basename "${out_dir}") (${model_path})"
  mkdir -p "${out_dir}"
  local out_abs
  out_abs="$(cd "${out_dir}" && pwd)"

  docker run --rm --platform linux/amd64 \
    -v "${out_abs}:/out" \
    -e MODEL_PATH="${model_path}" \
    "${docker_img}" \
    bash -c '
      pip install --quiet --no-warn-script-location keras 2>&1 | tail -1
      python3 -c "
import keras, numpy as np, os
m = keras.models.load_model(os.environ[\"MODEL_PATH\"], compile=False)
for i, layer in enumerate(m.layers):
    if not layer.weights:
        continue
    for w in layer.weights:
        arr = w.numpy()
        kind = \"kernel\" if \"kernel\" in w.name else (\"bias\" if \"bias\" in w.name else w.name)
        np.save(f\"/out/layer_{i}_{kind}.npy\", arr)
        print(f\"  layer_{i}_{kind}.npy  shape={arr.shape}  dtype={arr.dtype}\")
"
'
  log "done → ${out_dir}"
  ls -la "${out_dir}" | awk 'NR>1{print "    " $NF " (" $5 " B)"}'
}

# try_extract_small_from_image <out_dir> <model_path> <docker_image>
# Probes whether <model_path> exists in <docker_image> before extracting.
# Silently skips if absent (some trio/somatic small models are optional).
try_extract_small_from_image() {
  local out_dir="$1"
  local model_path="$2"
  local docker_img="$3"

  if [[ -d "${out_dir}" && -n "$(ls -A "${out_dir}" 2>/dev/null)" ]]; then
    log "SKIP (exists) $(basename "${out_dir}")"
    return 0
  fi

  log "Probing ${model_path} in ${docker_img##*/}..."
  if docker run --rm --platform linux/amd64 \
      "${docker_img}" \
      test -f "${model_path}" 2>/dev/null; then
    extract_small_from_image "${out_dir}" "${model_path}" "${docker_img}"
  else
    log "SKIP (not present in Docker) ${model_path}"
  fi
}

# ── Pull Docker images (once, idempotent) ─────────────────────────────────
log "Pulling Docker images (linux/amd64)..."
docker pull --platform linux/amd64 "${IMG_DV}"
docker pull --platform linux/amd64 "${IMG_DT}"
docker pull --platform linux/amd64 "${IMG_DS}"

# ══════════════════════════════════════════════════════════════════════════
# 1. DeepVariant germline big models — google/deepvariant:1.10.0
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepVariant germline big models ---"

# wgs — already present; included for completeness / re-extraction safety.
extract_dvw \
  "${MODELS_DIR}/wgs" \
  "${DVW_OUT}/wgs.dvw"

# wes
extract_dvw \
  "${MODELS_DIR}/wes" \
  "${DVW_OUT}/wes.dvw"

# pacbio
extract_dvw \
  "${MODELS_DIR}/pacbio" \
  "${DVW_OUT}/pacbio.dvw"

# ont
extract_dvw \
  "${MODELS_DIR}/ont" \
  "${DVW_OUT}/ont.dvw"

# hybrid (HYBRID_PACBIO_ILLUMINA)
extract_dvw \
  "${MODELS_DIR}/hybrid" \
  "${DVW_OUT}/hybrid.dvw"

# masseq
extract_dvw \
  "${MODELS_DIR}/masseq" \
  "${DVW_OUT}/masseq.dvw"

# rnaseq
extract_dvw \
  "${MODELS_DIR}/rnaseq" \
  "${DVW_OUT}/rnaseq.dvw"

# ══════════════════════════════════════════════════════════════════════════
# 2. DeepVariant germline small models — google/deepvariant:1.10.0
# Small model paths inside Docker: /opt/smallmodels/<variant>/model.keras
# WES, Hybrid, MaSeq, RNASeq have no trained_small_model_path.
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepVariant germline small models ---"

# wgs small model — already present; included for safety.
extract_small \
  "${DVW_OUT}/wgs_small_weights" \
  "/opt/smallmodels/wgs/model.keras" \
  "${DV_VERSION}"

# pacbio small model
extract_small \
  "${DVW_OUT}/pacbio_small_weights" \
  "/opt/smallmodels/pacbio/model.keras" \
  "${DV_VERSION}"

# ont small model
extract_small \
  "${DVW_OUT}/ont_small_weights" \
  "/opt/smallmodels/ont/model.keras" \
  "${DV_VERSION}"

# ══════════════════════════════════════════════════════════════════════════
# 3. DeepTrio big models — google/deepvariant:deeptrio-1.10.0
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepTrio big models ---"

# wgs child + parent — already present; included for safety.
extract_dvw \
  "${MODELS_DIR}/deeptrio.wgs_child" \
  "${DVW_OUT}/deeptrio.wgs_child.dvw" \
  "${IMG_DT}"

extract_dvw \
  "${MODELS_DIR}/deeptrio.wgs_parent" \
  "${DVW_OUT}/deeptrio.wgs_parent.dvw" \
  "${IMG_DT}"

# wes child + parent
extract_dvw \
  "${MODELS_DIR}/deeptrio.wes_child" \
  "${DVW_OUT}/deeptrio.wes_child.dvw" \
  "${IMG_DT}"

extract_dvw \
  "${MODELS_DIR}/deeptrio.wes_parent" \
  "${DVW_OUT}/deeptrio.wes_parent.dvw" \
  "${IMG_DT}"

# pacbio child + parent
extract_dvw \
  "${MODELS_DIR}/deeptrio.pacbio_child" \
  "${DVW_OUT}/deeptrio.pacbio_child.dvw" \
  "${IMG_DT}"

extract_dvw \
  "${MODELS_DIR}/deeptrio.pacbio_parent" \
  "${DVW_OUT}/deeptrio.pacbio_parent.dvw" \
  "${IMG_DT}"

# ont child + parent
extract_dvw \
  "${MODELS_DIR}/deeptrio.ont_child" \
  "${DVW_OUT}/deeptrio.ont_child.dvw" \
  "${IMG_DT}"

extract_dvw \
  "${MODELS_DIR}/deeptrio.ont_parent" \
  "${DVW_OUT}/deeptrio.ont_parent.dvw" \
  "${IMG_DT}"

# ══════════════════════════════════════════════════════════════════════════
# 4. DeepTrio small models — google/deepvariant:deeptrio-1.10.0
# WGS paths confirmed in extract_small_model_weights.sh comments:
#   /opt/smallmodels/deeptrio/wgs/child/model.keras
#   /opt/smallmodels/deeptrio/wgs/parent/model.keras
# WES/PacBio/ONT paths probed at runtime (may not exist upstream).
# Note: extract_small_model_weights.sh uses google/deepvariant:<tag>, so
# "deeptrio-1.10.0" maps to google/deepvariant:deeptrio-1.10.0 — correct.
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepTrio small models ---"

# wgs child + parent — already present; included for safety.
extract_small \
  "${DVW_OUT}/deeptrio_wgs_child_small" \
  "/opt/smallmodels/deeptrio/wgs/child/model.keras" \
  "deeptrio-${DV_VERSION}"

extract_small \
  "${DVW_OUT}/deeptrio_wgs_parent_small" \
  "/opt/smallmodels/deeptrio/wgs/parent/model.keras" \
  "deeptrio-${DV_VERSION}"

# wes/pacbio/ont: probe first since upstream may not ship these.
try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_wes_child_small" \
  "/opt/smallmodels/deeptrio/wes/child/model.keras" \
  "${IMG_DT}"

try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_wes_parent_small" \
  "/opt/smallmodels/deeptrio/wes/parent/model.keras" \
  "${IMG_DT}"

try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_pacbio_child_small" \
  "/opt/smallmodels/deeptrio/pacbio/child/model.keras" \
  "${IMG_DT}"

try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_pacbio_parent_small" \
  "/opt/smallmodels/deeptrio/pacbio/parent/model.keras" \
  "${IMG_DT}"

try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_ont_child_small" \
  "/opt/smallmodels/deeptrio/ont/child/model.keras" \
  "${IMG_DT}"

try_extract_small_from_image \
  "${DVW_OUT}/deeptrio_ont_parent_small" \
  "/opt/smallmodels/deeptrio/ont/parent/model.keras" \
  "${IMG_DT}"

# ══════════════════════════════════════════════════════════════════════════
# 5. DeepSomatic big models — google/deepsomatic:1.10.0
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepSomatic big models ---"

# wgs — already present; included for safety.
extract_dvw \
  "${MODELS_DIR}/deepsomatic.wgs" \
  "${DVW_OUT}/deepsomatic.wgs.dvw" \
  "${IMG_DS}"

# wes
extract_dvw \
  "${MODELS_DIR}/deepsomatic.wes" \
  "${DVW_OUT}/deepsomatic.wes.dvw" \
  "${IMG_DS}"

# pacbio
extract_dvw \
  "${MODELS_DIR}/deepsomatic.pacbio" \
  "${DVW_OUT}/deepsomatic.pacbio.dvw" \
  "${IMG_DS}"

# ont
extract_dvw \
  "${MODELS_DIR}/deepsomatic.ont" \
  "${DVW_OUT}/deepsomatic.ont.dvw" \
  "${IMG_DS}"

# ffpe_wgs
extract_dvw \
  "${MODELS_DIR}/deepsomatic.ffpe_wgs" \
  "${DVW_OUT}/deepsomatic.ffpe_wgs.dvw" \
  "${IMG_DS}"

# ffpe_wes
extract_dvw \
  "${MODELS_DIR}/deepsomatic.ffpe_wes" \
  "${DVW_OUT}/deepsomatic.ffpe_wes.dvw" \
  "${IMG_DS}"

# ══════════════════════════════════════════════════════════════════════════
# 5b. DeepSomatic tumor-only big models — google/deepsomatic:1.10.0
# Tumor-only models use separate SavedModels: 8 channels for WGS/WES/FFPE
# (base-6 + insert_size + allele_frequency), 10 for PacBio/ONT (base-6 +
# haplotype + suppl/fuzzy + alt_aligned×2 + allele_frequency), h=100.
# No small models exist for any tumor-only variant.
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepSomatic tumor-only big models ---"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.wgs_tumor_only" \
  "${DVW_OUT}/deepsomatic.wgs_tumor_only.dvw" \
  "${IMG_DS}"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.wes_tumor_only" \
  "${DVW_OUT}/deepsomatic.wes_tumor_only.dvw" \
  "${IMG_DS}"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.ffpe_wgs_tumor_only" \
  "${DVW_OUT}/deepsomatic.ffpe_wgs_tumor_only.dvw" \
  "${IMG_DS}"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.ffpe_wes_tumor_only" \
  "${DVW_OUT}/deepsomatic.ffpe_wes_tumor_only.dvw" \
  "${IMG_DS}"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.pacbio_tumor_only" \
  "${DVW_OUT}/deepsomatic.pacbio_tumor_only.dvw" \
  "${IMG_DS}"

extract_dvw \
  "${MODELS_DIR}/deepsomatic.ont_tumor_only" \
  "${DVW_OUT}/deepsomatic.ont_tumor_only.dvw" \
  "${IMG_DS}"

# ══════════════════════════════════════════════════════════════════════════
# 5c. DeepSomatic Panel-of-Normals (PON) — google/deepsomatic:1.10.0
# population_vcfs in tumor-only example_info.json flags_for_calling.
# Required for allele-frequency filtering in make_examples tumor-only mode.
# Extracted once; stored alongside model weights.
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepSomatic PON files (Illumina + PacBio/ONT) ---"

PON_OUT="${DVW_OUT}/deepsomatic_pon"
mkdir -p "${PON_OUT}"
PON_ABS="$(cd "${PON_OUT}" && pwd)"

# Illumina PON — used by WGS/WES/FFPE tumor-only modes
PON_ILMN="${PON_OUT}/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz"
if [[ -s "${PON_ILMN}" ]]; then
  log "SKIP (exists) AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz"
else
  log "Extracting Illumina PON from ${IMG_DS}"
  docker run --rm --platform linux/amd64 \
    -v "${PON_ABS}:/out" \
    "${IMG_DS}" \
    bash -c '
      src=/opt/models/deepsomatic/pons/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz
      if [[ -f "${src}" ]]; then
        cp "${src}" /out/
        [[ -f "${src}.tbi" ]] && cp "${src}.tbi" /out/ || \
          tabix -p vcf "/out/$(basename ${src})" 2>/dev/null || true
        echo "done"
      else
        echo "WARNING: Illumina PON not found at ${src}" >&2
      fi
    '
fi

# PacBio/ONT PON — used by PacBio and ONT tumor-only modes (CoLoRSdb ~254 MB)
PON_PB="${PON_OUT}/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz"
if [[ -s "${PON_PB}" ]]; then
  log "SKIP (exists) AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz"
else
  log "Extracting PacBio/ONT PON from ${IMG_DS}"
  docker run --rm --platform linux/amd64 \
    -v "${PON_ABS}:/out" \
    "${IMG_DS}" \
    bash -c '
      src=/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz
      if [[ -f "${src}" ]]; then
        cp "${src}" /out/
        [[ -f "${src}.tbi" ]] && cp "${src}.tbi" /out/ || \
          tabix -p vcf "/out/$(basename ${src})" 2>/dev/null || true
        echo "done"
      else
        echo "WARNING: PacBio PON not found at ${src}" >&2
      fi
    '
fi

# ══════════════════════════════════════════════════════════════════════════
# 6. DeepSomatic small models — google/deepsomatic:1.10.0
# Paths from example_info.json trained_small_model_path:
#   deepsomatic.wgs     → /opt/smallmodels/wgs/model.keras
#   deepsomatic.pacbio  → /opt/smallmodels/pacbio/model.keras
#   deepsomatic.ont     → /opt/smallmodels/ont/model.keras
#   deepsomatic.ffpe_wgs → /opt/smallmodels/ffpe_wgs/model.keras
#   deepsomatic.wes     → no small model
#   deepsomatic.ffpe_wes → no small model
# Uses IMG_DS (google/deepsomatic:1.10.0) directly, not extract_small_model_weights.sh
# which is hardcoded to google/deepvariant:<tag>.
# ══════════════════════════════════════════════════════════════════════════

log "--- DeepSomatic small models ---"

try_extract_small_from_image \
  "${DVW_OUT}/deepsomatic_wgs_small" \
  "/opt/smallmodels/wgs/model.keras" \
  "${IMG_DS}"

try_extract_small_from_image \
  "${DVW_OUT}/deepsomatic_pacbio_small" \
  "/opt/smallmodels/pacbio/model.keras" \
  "${IMG_DS}"

try_extract_small_from_image \
  "${DVW_OUT}/deepsomatic_ont_small" \
  "/opt/smallmodels/ont/model.keras" \
  "${IMG_DS}"

try_extract_small_from_image \
  "${DVW_OUT}/deepsomatic_ffpe_wgs_small" \
  "/opt/smallmodels/ffpe_wgs/model.keras" \
  "${IMG_DS}"

# ══════════════════════════════════════════════════════════════════════════
# 7. Pangenome-aware DeepVariant — already present; included for safety.
#    No small model (no trained_small_model_path in example_info.json).
# ══════════════════════════════════════════════════════════════════════════

log "--- Pangenome-aware DeepVariant ---"

extract_dvw \
  "${MODELS_DIR}/pangenome_aware_deepvariant" \
  "${DVW_OUT}/pangenome.wgs.dvw"

# ══════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════

log "All extractions complete. Output directory: ${DVW_OUT}"
echo ""
echo "  .dvw files:"
ls -lh "${DVW_OUT}"/*.dvw 2>/dev/null \
  | awk '{printf "    %-48s %s\n", $NF, $5}' \
  || echo "    (none)"
echo ""
echo "  small-model weight directories:"
for d in "${DVW_OUT}"/*_small; do
  [[ -d "${d}" ]] || continue
  count=$(ls "${d}"/*.npy 2>/dev/null | wc -l | tr -d ' ')
  printf "    %-48s %s .npy files\n" "$(basename "${d}")" "${count}"
done
