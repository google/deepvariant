#!/usr/bin/env bash
# Fetch + convert every DeepVariant model variant (big + small) to a Mac
# arm64 .mlpackage, ready to be shipped by the deepvariant-models Homebrew
# formula.
#
# Variants:
#   DeepVariant (germline, single-sample):
#     wgs, wes, pacbio, ont_r104, hybrid_pacbio_illumina, masseq, rnaseq
#   DeepTrio (germline, mother/father/child):
#     deeptrio.wgs.{child,parent}, deeptrio.pacbio.{child,parent},
#     deeptrio.ont.{child,parent}, deeptrio.wes.{child,parent}
#   DeepSomatic (tumor/normal):
#     deepsomatic.{wgs,wes,ffpe_wgs,ffpe_wes,ont,pacbio}
#   Pangenome-aware DeepVariant (12-channel input):
#     pangenome_wgs
#
# Each variant ships TWO models: a big InceptionV3 (image pileup) and a
# small MLP (70 features). The small model lives at
# /opt/smallmodels/<variant>/model.keras inside google/deepvariant:1.10.0;
# the big one lives in gs://deepvariant/models/DeepVariant/1.10.0/savedmodels/.
#
# Output: tools/conversion/models/<name>.mlpackage and
#         tools/conversion/models/<name>_small.mlpackage

set -euo pipefail
cd "$(dirname "$0")"

DV_VERSION="${DV_VERSION:-1.10.0}"

# ── DeepVariant ─────────────────────────────────────────────────────────────
DV_VARIANTS=(
  wgs
  wes
  pacbio
  ont_r104
  hybrid_pacbio_illumina
)

# ── DeepTrio (3 sub-models per variant: parent, child) ──────────────────────
# Note: trio uses 6-channel input (different shape from DV's 7-channel).
DEEPTRIO_VARIANTS=(
  # disabled by default — uncomment to convert.
  # deeptrio.wgs.parent
  # deeptrio.wgs.child
  # deeptrio.pacbio.parent
  # deeptrio.pacbio.child
)

# ── DeepSomatic ─────────────────────────────────────────────────────────────
DEEPSOMATIC_VARIANTS=(
  # deepsomatic.wgs
  # deepsomatic.wes
)

convert_one() {
  local NAME="$1"
  ./fetch_savedmodel.sh "${NAME}" || return 1
  ./convert_via_docker.sh "models/${NAME}" "models/${NAME}.mlpackage"
  if docker run --rm --platform linux/amd64 \
      "google/deepvariant:${DV_VERSION}" \
      test -f "/opt/smallmodels/${NAME}/model.keras" 2>/dev/null; then
    ./convert_small_model.sh "${NAME}"
  else
    echo "==> no /opt/smallmodels/${NAME}/model.keras (skipping small model)"
  fi
}

for v in "${DV_VARIANTS[@]}"; do
  echo "============================================================"
  echo "== Converting ${v}"
  echo "============================================================"
  convert_one "${v}" || echo "  (conversion of ${v} failed; continuing)"
done

for v in "${DEEPTRIO_VARIANTS[@]}"; do
  echo "============================================================"
  echo "== DeepTrio: ${v}"
  echo "============================================================"
  echo "  NOTE: Trio uses 6-channel input. The conversion script's input"
  echo "        shape (100,221,7) needs to be adjusted to (100,221,6) for"
  echo "        trio big-model conversion. See PORT_LOG for full plan."
  # convert_one "${v}" || echo "  (failed; continuing)"
done

for v in "${DEEPSOMATIC_VARIANTS[@]}"; do
  echo "============================================================"
  echo "== DeepSomatic: ${v}"
  echo "============================================================"
  # convert_one "${v}" || echo "  (failed; continuing)"
done

echo
echo "==> Conversion complete. Models in tools/conversion/models/:"
ls -d models/*.mlpackage 2>/dev/null || echo "  (none yet)"
