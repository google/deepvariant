#!/usr/bin/env bash
# Extract small_model weights from an upstream Docker image into a
# directory of FP32 NumPy `.npy` files. The runtime
# `SmallModel::Load(<dir>)` reads them and runs a deterministic
# BNNS-CPU FP32 MLP (Phase 5.5d/7) — bit-equal to Docker's TF/Keras
# inference.
#
# Output: <out_dir>/{layer_0_kernel,layer_0_bias,layer_1_kernel,
#                    layer_1_bias,layer_2_kernel,layer_2_bias}.npy
# (~2.4 MB total for the standard 3-layer MLP).
#
# Usage:
#   tools/conversion/extract_small_model_weights.sh <out_dir> [model_path] [docker_tag]
#
# Examples (Phase 5.5/6 model variants):
#   tools/conversion/extract_small_model_weights.sh validation/work/wgs_small_weights
#       (defaults to /opt/smallmodels/wgs/model.keras in google/deepvariant:1.10.0)
#
#   tools/conversion/extract_small_model_weights.sh validation/work/deeptrio_wgs_child_small \
#       /opt/smallmodels/deeptrio/wgs/child/model.keras  deeptrio-1.10.0
#
#   tools/conversion/extract_small_model_weights.sh validation/work/deeptrio_wgs_parent_small \
#       /opt/smallmodels/deeptrio/wgs/parent/model.keras  deeptrio-1.10.0
#
# Re-run whenever the upstream Docker image bumps weights. Stable
# across DV 1.10.0 sub-releases.

set -euo pipefail

OUT_DIR="${1:?usage: $0 <out_dir> [model_path] [docker_tag]}"
MODEL_PATH="${2:-/opt/smallmodels/wgs/model.keras}"
DV_TAG="${3:-1.10.0}"

mkdir -p "${OUT_DIR}"
OUT_ABS="$(cd "${OUT_DIR}" && pwd)"

echo "==> Extracting ${MODEL_PATH} from google/deepvariant:${DV_TAG}"
echo "    out: ${OUT_ABS}"

docker run --rm --platform linux/amd64 \
  -v "${OUT_ABS}:/out" \
  -e MODEL_PATH="${MODEL_PATH}" \
  "google/deepvariant:${DV_TAG}" \
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

echo "==> done"
ls -la "${OUT_ABS}" | awk 'NR>1{print "    " $NF " (" $5 " B)"}'
