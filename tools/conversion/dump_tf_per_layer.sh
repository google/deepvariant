#!/usr/bin/env bash
# Dump per-layer reference outputs for a DeepVariant SavedModel via the
# upstream Docker image's TF 2.16. Produces `<out_dir>/<tap>.npy` for
# every Metal tap (`stem_s1a`, `stem_s2a`, …, `7c`, `gap`) plus a
# `_savedmodel_softmax.npy` (canonical Docker output) and a
# `_reimpl_softmax.npy` (TF reimpl sanity check).
#
# Mirrors convert_via_docker.sh's pattern (linux/amd64 emulated, RO
# input mount, RW output mount, tools/conversion mounted at /work so
# the script can `import` shared helpers).
#
# Usage:
#   ./dump_tf_per_layer.sh <savedmodel_dir> <out_dir>             # seed-0
#   ./dump_tf_per_layer.sh <savedmodel_dir> <out_dir> <input.npy> # real input

set -euo pipefail

MODEL_DIR="${1:?usage: $0 <savedmodel_dir> <out_dir> [input.npy]}"
OUT_DIR="${2:?usage: $0 <savedmodel_dir> <out_dir> [input.npy]}"
INPUT_NPY="${3:-}"
DV_VERSION="${DV_VERSION:-1.10.0}"

if [[ ! -f "${MODEL_DIR}/saved_model.pb" ]]; then
  echo "error: ${MODEL_DIR} does not contain saved_model.pb" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
MODEL_DIR_ABS="$(cd "${MODEL_DIR}" && pwd)"
OUT_DIR_ABS="$(cd "${OUT_DIR}" && pwd)"
WORK_DIR_ABS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==> Dumping per-layer references"
echo "    model:  ${MODEL_DIR_ABS}"
echo "    out:    ${OUT_DIR_ABS}"
echo "    image:  google/deepvariant:${DV_VERSION}"

if [[ -n "${INPUT_NPY}" ]]; then
  if [[ ! -f "${INPUT_NPY}" ]]; then
    echo "error: input.npy not found: ${INPUT_NPY}" >&2
    exit 1
  fi
  INPUT_DIR_ABS="$(cd "$(dirname "${INPUT_NPY}")" && pwd)"
  INPUT_NAME="$(basename "${INPUT_NPY}")"
  echo "    input:  ${INPUT_DIR_ABS}/${INPUT_NAME}"
  docker run --rm --platform linux/amd64 \
    -v "${MODEL_DIR_ABS}:/in:ro" \
    -v "${OUT_DIR_ABS}:/out" \
    -v "${WORK_DIR_ABS}:/work:ro" \
    -v "${INPUT_DIR_ABS}:/inp:ro" \
    "google/deepvariant:${DV_VERSION}" \
    python3 /work/dump_tf_per_layer.py /in /out "/inp/${INPUT_NAME}"
else
  docker run --rm --platform linux/amd64 \
    -v "${MODEL_DIR_ABS}:/in:ro" \
    -v "${OUT_DIR_ABS}:/out" \
    -v "${WORK_DIR_ABS}:/work:ro" \
    "google/deepvariant:${DV_VERSION}" \
    python3 /work/dump_tf_per_layer.py /in /out
fi

echo "==> done"
