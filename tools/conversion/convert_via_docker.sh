#!/usr/bin/env bash
# Convert a DeepVariant SavedModel → Core ML .mlpackage via the upstream
# Docker image's coremltools. Shape is auto-detected from the SavedModel's
# model.example_info.json so the same recipe works for WGS (7 ch),
# PacBio/ONT (10 ch), trio (6 ch), pangenome (12 ch), etc.
#
# Bit-parity verified vs upstream call_variants on a chr20 fixture for
# WGS: 100.000% argmax + softmax max-abs = 0.000000 (gate ≤1e-3).
#
# Usage: ./convert_via_docker.sh <model_dir> <out.mlpackage>

set -euo pipefail

MODEL_DIR="${1:?usage: $0 <savedmodel_dir> <out.mlpackage>}"
OUT="${2:?usage: $0 <savedmodel_dir> <out.mlpackage>}"
DV_VERSION="${DV_VERSION:-1.10.0}"

if [[ ! -f "${MODEL_DIR}/saved_model.pb" ]]; then
  echo "error: ${MODEL_DIR} does not contain saved_model.pb" >&2
  exit 1
fi

# Auto-detect shape from model.example_info.json. Falls back to WGS shape.
if [[ -f "${MODEL_DIR}/model.example_info.json" ]]; then
  SHAPE_JSON=$(python3 -c "
import json, sys
d = json.load(open('${MODEL_DIR}/model.example_info.json'))
print(','.join(str(x) for x in d['shape']))
")
  [[ -n "${SHAPE_JSON}" && "${SHAPE_JSON}" == *,*,* ]] || { echo "error: failed to parse shape from ${MODEL_DIR}/model.example_info.json" >&2; exit 1; }
else
  SHAPE_JSON="100,221,7"
fi
H=$(echo "$SHAPE_JSON" | cut -d, -f1)
W=$(echo "$SHAPE_JSON" | cut -d, -f2)
C=$(echo "$SHAPE_JSON" | cut -d, -f3)

MODEL_DIR_ABS="$(cd "${MODEL_DIR}" && pwd)"
OUT_DIR="$(cd "$(dirname "${OUT}")" && pwd)"
OUT_NAME="$(basename "${OUT}")"

echo "==> Converting ${MODEL_DIR_ABS} → ${OUT_DIR}/${OUT_NAME}"
echo "    Auto-detected input shape: (N, ${H}, ${W}, ${C})"

docker run --rm --platform linux/amd64 \
  -v "${MODEL_DIR_ABS}:/in:ro" \
  -v "${OUT_DIR}:/out" \
  -e "H=${H}" -e "W=${W}" -e "C=${C}" -e "OUT_NAME=${OUT_NAME}" \
  "google/deepvariant:${DV_VERSION}" \
  bash -c "pip install --quiet coremltools==7.2 2>&1 | tail -1 && \
    python3 -c '
import os, coremltools as ct, numpy as np
H = int(os.environ[\"H\"])
W = int(os.environ[\"W\"])
C = int(os.environ[\"C\"])
out_name = os.environ[\"OUT_NAME\"]
print(f\"Converting at FLOAT32 with shape (N, {H}, {W}, {C})...\")
mlmodel = ct.convert(
    \"/in\",
    convert_to=\"mlprogram\",
    source=\"tensorflow\",
    inputs=[ct.TensorType(
        name=\"input_1\",
        shape=(ct.RangeDim(1, 4096), H, W, C),
        dtype=np.float32,
    )],
    outputs=[ct.TensorType(name=\"Identity\", dtype=np.float32)],
    minimum_deployment_target=ct.target.macOS14,
    compute_precision=ct.precision.FLOAT32,
)
mlmodel.save(f\"/out/{out_name}\")
print(f\"Saved /out/{out_name}\")
'"

echo "==> done"
