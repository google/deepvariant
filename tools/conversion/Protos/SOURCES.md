# Vendored protobuf sources

All `.proto` files under `Protos/tensorflow/` are vendored verbatim from upstream and **not patched**. Their license is each upstream project's own. Re-fetch via the commands below if anything changes upstream.

## TensorFlow — `tensorflow/r2.16` branch (Apache-2.0)

Source: <https://github.com/tensorflow/tensorflow>
Branch: `r2.16` (matches the TF version that DeepVariant 1.10 SavedModels were written by).
Fetched: 2026-04-25.

Files (26 — `error_codes.proto` and `debug_event.proto` were dropped because they pull in `tsl/protobuf/error_codes.proto` from a separate Google package, and neither is needed for SavedModel parsing):

```text
core/framework/allocation_description.proto
core/framework/attr_value.proto
core/framework/cost_graph.proto
core/framework/device_attributes.proto
core/framework/full_type.proto
core/framework/function.proto
core/framework/graph.proto
core/framework/graph_debug_info.proto
core/framework/node_def.proto
core/framework/op_def.proto
core/framework/resource_handle.proto
core/framework/step_stats.proto
core/framework/tensor.proto
core/framework/tensor_description.proto
core/framework/tensor_shape.proto
core/framework/tensor_slice.proto
core/framework/types.proto
core/framework/variable.proto
core/framework/versions.proto
core/protobuf/meta_graph.proto
core/protobuf/saved_model.proto
core/protobuf/saved_object_graph.proto
core/protobuf/saver.proto
core/protobuf/struct.proto
core/protobuf/tensor_bundle.proto
core/protobuf/trackable_object_graph.proto
```

Re-fetch:

```sh
TF_REF="r2.16"
BASE="https://raw.githubusercontent.com/tensorflow/tensorflow/${TF_REF}/tensorflow"
cd tools/conversion/Protos/tensorflow
for f in <list above>; do
  curl -fsSL -o "${f}" "${BASE}/${f}"
done
```

## Generation

Python bindings are generated under `tools/conversion/Generated/` (gitignored):

```sh
cd tools/conversion
rm -rf Generated && mkdir Generated
protoc --python_out=Generated/ \
       --proto_path=Protos \
       $(cd Protos && find tensorflow -name '*.proto')
touch Generated/__init__.py
```

After generation, the `tensorflow.core.protobuf.*_pb2` and `tensorflow.core.framework.*_pb2` modules become importable when `Generated/` is on the Python path:

```python
import sys; sys.path.insert(0, "tools/conversion/Generated")
from tensorflow.core.protobuf import saved_model_pb2, tensor_bundle_pb2
```

The proto_path must be the parent of `tensorflow/`, not `tensorflow/` itself, because the proto files use absolute-style imports like `import "tensorflow/core/framework/graph.proto"`.

Bindings are regenerated on demand by `setup_venvs.sh` (or the snippet above).

## Why we vendor instead of pip-install

The natural way to get TF's `.proto` definitions is `pip install tensorflow`, which we explicitly forbid (Voie B refined — TF banned in v2). Vendoring the schema files alone is ~110 KB and gives us proto bindings via `protoc --python_out` with no TF runtime.
