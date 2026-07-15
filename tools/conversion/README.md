# Phase 0 — Inference framework bench (dev-time Python tooling)

> **Dev-time only — never shipped to users.** Python lives here because `coremltools` and `mlx` are Python-only Apple packages, and coremltools' MIL builder carries years of edge-case maturity (BatchNorm fused vs not, ConcatV2 axis handling, NHWC↔NCHW layout, FusedBatchNormV3 epsilon, etc.) that we don't want to re-implement in Swift. The conversion pipeline is **TF-free**: weights are read straight from the SavedModel's TensorBundle by the pure-protobuf readers in this directory (`savedmodel_reader.py`, `tensor_bundle_reader.py`), with no TensorFlow, PyTorch, ONNX, or tensorflow-metal involved. The user-facing `deepvariant` binary stays 100 % C++/Obj-C++ with no embedded Python interpreter (verified by `otool -L` in Phase 5). See `CLAUDE.md` for the dev-time vs runtime split.

Two-way A/B bench (Core ML vs MLX) on the real WGS SavedModel: produces latency, throughput, GPU/ANE residency, and parity-vs-Linux measurements that feed the Phase 0 ADR. (The earlier tensorflow-metal voie was dropped — tensorflow-metal 1.2.0 has been unmaintained since mid-2024.)

## Layout

```text
tools/conversion/
├── .python-version            # pyenv pin: 3.11.x
├── requirements-coreml.txt    # coremltools 9.0+ (TF-free MIL path)
├── requirements-mlx.txt       # MLX 0.21+ + safetensors (TF-free)
├── setup_venvs.sh             # creates venv-coreml / venv-mlx (TF-free)
├── fetch_savedmodel.sh        # pulls gs://deepvariant/models/DeepVariant/1.10.0/<name>
├── savedmodel_reader.py       # pure-protobuf SavedModel reader (no TF)
├── tensor_bundle_reader.py    # pure-protobuf TensorBundle weight reader (no TF)
├── convert_coreml.py          # TensorBundle weights -> .mlpackage (coremltools MIL)
├── convert_mlx.py             # SavedModel weights -> safetensors (MLX-friendly)
├── bench.py                   # latency + powermetrics GPU residency + softmax capture
├── parity_check.py            # softmax max-abs / argmax disagreement vs reference
└── models/                    # gitignored, where SavedModels and outputs live
```

## Pipeline

```sh
# one-time setup (pyenv install 3.11 + two venv pip installs)
./setup_venvs.sh

# pull the real WGS SavedModel (~700 MB)
./fetch_savedmodel.sh wgs

# convert each voie
source venv-coreml/bin/activate
python convert_coreml.py --bundle models/wgs/variables/variables --output models/wgs.mlpackage
deactivate

source venv-mlx/bin/activate
python convert_mlx.py --saved-model models/wgs --output models/wgs.mlx.safetensors
deactivate

# bench each on the same 1000-example reference set
source venv-coreml/bin/activate
python bench.py --backend coreml --model models/wgs.mlpackage \
  --examples ../reference/cache/wgs_chr20_1000.tfrecord \
  --output ../../benchmarks/coreml_wgs.json \
  --output-cv ../../benchmarks/coreml_wgs.cv.tfrecord
deactivate
# (repeat for mlx)

# parity vs Linux reference
python parity_check.py \
  --reference ../reference/output/wgs/call_variants_chr20.tfrecord \
  --candidates ../../benchmarks/coreml_wgs.cv.tfrecord \
              ../../benchmarks/mlx_wgs.cv.tfrecord
```

## Why two pinned venvs

Core ML and MLX have incompatible Python dependency requirements, so each gets its own venv. Both are **TF-free** — `setup_venvs.sh` hard-fails if `tensorflow` is importable in either.

| venv | Python | Key deps |
| --- | --- | --- |
| venv-coreml | 3.11 | coremltools 9.0+ (direct MIL path; no TF SavedModel conversion) |
| venv-mlx | 3.11 | MLX 0.21+, safetensors (Inception-v3 rebuilt in MLX at bench time) |

Both venvs pin `numpy < 2` (see the requirements files).

## Stop conditions

- If the parity check shows argmax disagreement on **any** of the 1000 examples, the framework is rejected — no exceptions.
- If `bench.py` reports `gpu_power=0` AND `ane_power=0` for a backend, that backend has a config bug, not a perf result.
