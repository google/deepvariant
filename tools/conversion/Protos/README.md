# Vendored TF `.proto` files (TF-free SavedModel reading)

26 `.proto` files vendored from `tensorflow/r2.16` (Apache-2.0). See `SOURCES.md` for the exact list, the upstream branch, and the re-fetch commands.

The vendored set is just the schema definitions — no TensorFlow runtime, no Python TF package. We compile them with system `protoc --python_out=Generated/` and import the resulting `*_pb2` modules from `savedmodel_reader.py`.

## Generate Python bindings

```sh
cd tools/conversion
mkdir -p Generated
protoc --python_out=Generated/ \
       --proto_path=Protos/tensorflow \
       $(find Protos/tensorflow -name '*.proto')
touch Generated/__init__.py
```

After generation, all `core/protobuf/*.proto` are accessible as `core.protobuf.*_pb2`, and `core/framework/*.proto` as `core.framework.*_pb2`.

The `Generated/` directory is `.gitignore`d — bindings are regenerated on every `setup_venvs.sh` run.
