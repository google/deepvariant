#!/usr/bin/env bash
# Convert the upstream "small model" (a 70-feature → 3-class MLP) from
# Keras .keras format to a Core ML .mlpackage.
#
# The small model is what produces upstream's tight QUAL/GQ on ~84 % of
# WGS candidates: a 3-layer dense network (70 → 750 → 750 → 3) running
# on hand-engineered features (allele counts, qualities, VAF context).
# The big InceptionV3 model only kicks in for the remaining ~16 % that
# don't pass the small-model GQ threshold (snp=20, indel=28).
#
# Usage: ./convert_small_model.sh <variant>
# e.g.   ./convert_small_model.sh wgs       (→ models/wgs_small.mlpackage)
# Variants ship in /opt/smallmodels/<variant>/model.keras inside the
# google/deepvariant:1.10.0 Docker image.

set -euo pipefail
VARIANT="${1:?usage: $0 <wgs|wes|pacbio|ont_r104>}"
DV_VERSION="${DV_VERSION:-1.10.0}"

OUT_DIR="$(cd "$(dirname "$0")" && pwd)/models"
mkdir -p "${OUT_DIR}/${VARIANT}_small"

echo "==> Extracting upstream /opt/smallmodels/${VARIANT}/ ..."
docker run --rm --platform linux/amd64 \
  -v "${OUT_DIR}/${VARIANT}_small:/copy" \
  "google/deepvariant:${DV_VERSION}" \
  bash -c "cp -r /opt/smallmodels/${VARIANT}/* /copy/"

echo "==> Converting model.keras → ${VARIANT}_small.mlpackage ..."
docker run --rm --platform linux/amd64 \
  -v "${OUT_DIR}/${VARIANT}_small:/in" \
  -v "${OUT_DIR}:/out" \
  "google/deepvariant:${DV_VERSION}" \
  bash -c "
pip install --quiet keras==3.5.0 coremltools==7.2 2>&1 | tail -2
python3 -c '
import keras, tensorflow as tf, coremltools as ct, numpy as np

m = keras.models.load_model(\"/in/model.keras\", compile=False)
print(\"Loaded:\", m.count_params(), \"params\")

# Wrap as TF SavedModel for coremltools.
class Wrap(tf.Module):
    def __init__(self, k): self.k = k
    @tf.function(input_signature=[tf.TensorSpec((None, 70), tf.float32, name=\"input_1\")])
    def __call__(self, x): return {\"Identity\": self.k(x)}

w = Wrap(m)
tf.saved_model.save(w, \"/tmp/sm\", signatures=w.__call__.get_concrete_function())

mlmodel = ct.convert(
    \"/tmp/sm\",
    convert_to=\"mlprogram\",
    source=\"tensorflow\",
    inputs=[ct.TensorType(name=\"input_1\",
                          shape=(ct.RangeDim(1, 65536), 70), dtype=np.float32)],
    outputs=[ct.TensorType(name=\"Identity\", dtype=np.float32)],
    minimum_deployment_target=ct.target.macOS14,
    compute_precision=ct.precision.FLOAT32,
)
mlmodel.save(\"/out/${VARIANT}_small.mlpackage\")
'
"

echo "==> wrote ${OUT_DIR}/${VARIANT}_small.mlpackage"
