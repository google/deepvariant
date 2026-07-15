#!/usr/bin/env bash
# Set up the two Phase 0 conversion venvs (coreml / mlx). TF-free.
# tf-metal voie was dropped: tensorflow-metal 1.2.0 is unmaintained since mid-2024.
# Idempotent: safe to re-run.
set -euo pipefail

cd "$(dirname "$0")"

PYTHON_VERSION="$(cat .python-version)"

if ! command -v pyenv >/dev/null 2>&1; then
  echo "error: pyenv not found. brew install pyenv" >&2
  exit 1
fi

if ! pyenv versions --bare | grep -qx "${PYTHON_VERSION}"; then
  echo "==> installing Python ${PYTHON_VERSION} via pyenv"
  pyenv install "${PYTHON_VERSION}"
fi

PYBIN="$(pyenv root)/versions/${PYTHON_VERSION}/bin/python3"

setup_venv() {
  local name="$1"
  local req="requirements-${name}.txt"
  local venv="venv-${name}"

  if [[ -d "${venv}" ]]; then
    echo "==> ${venv} exists, refreshing pinned deps"
  else
    echo "==> creating ${venv}"
    "${PYBIN}" -m venv "${venv}"
  fi

  # shellcheck disable=SC1091
  source "${venv}/bin/activate"
  python -m pip install --upgrade pip wheel
  python -m pip install -r "${req}"
  python -c "import sys; print('python', sys.version)"

  # Hard guard: tensorflow MUST NOT be importable from any of our venvs.
  if python -c "import tensorflow" 2>/dev/null; then
    echo "FATAL: tensorflow is importable in ${venv} — TF-free policy violated" >&2
    deactivate
    exit 1
  fi

  case "${name}" in
    coreml)
      python -c "import torch, coremltools as ct; print('torch', torch.__version__, 'coremltools', ct.__version__)"
      ;;
    mlx)
      python -c "import mlx.core as mx; print('mlx default device:', mx.default_device())"
      ;;
  esac
  deactivate
}

setup_venv coreml
setup_venv mlx

echo "==> generating Python protobuf bindings under Generated/"
if ! command -v protoc >/dev/null 2>&1; then
  echo "error: protoc not found. brew install protobuf" >&2
  exit 1
fi
rm -rf Generated
mkdir -p Generated
( cd Protos && protoc --python_out=../Generated -I=. \
    $(find tensorflow -name '*.proto') )
touch Generated/__init__.py
echo "==> generated $(find Generated -name '*_pb2.py' | wc -l | tr -d ' ') Python modules"

echo "==> both venvs ready (TF-free)"
