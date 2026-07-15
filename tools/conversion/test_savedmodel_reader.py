"""Round-trip test: build a synthetic SavedModel proto in memory, write it,
parse it back via SavedModelReader. Run from `tools/conversion/`:

    python test_savedmodel_reader.py

No real SavedModel artefacts are needed — this exercises only the
proto-parsing half of the reader. The TensorBundle (weights) half is
implemented now; here we only assert that it raises on a malformed
(zeroed) variables.index rather than building a valid bundle inline.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# Bring the generated bindings onto sys.path for this dev test.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "Generated"))

from savedmodel_reader import SavedModelReader  # noqa: E402


def _build_synth_savedmodel(tmp: Path) -> Path:
    from tensorflow.core.protobuf import (  # type: ignore
        saved_model_pb2, meta_graph_pb2,
    )
    from tensorflow.core.framework import (  # type: ignore
        graph_pb2, node_def_pb2, attr_value_pb2,
        tensor_shape_pb2, types_pb2,
    )

    sm = saved_model_pb2.SavedModel(saved_model_schema_version=1)
    meta = sm.meta_graphs.add()

    # GraphDef with one VarHandleOp + one Conv2D + one Identity.
    g = meta.graph_def

    var = g.node.add(name="conv_kernel", op="VarHandleOp")
    var.attr["dtype"].type = types_pb2.DT_FLOAT
    var.attr["shared_name"].s = b"conv_kernel"

    conv = g.node.add(name="conv_out", op="Conv2D",
                      input=["input_image", "conv_kernel"])
    conv.attr["T"].type = types_pb2.DT_FLOAT
    conv.attr["padding"].s = b"SAME"
    conv.attr["strides"].list.i.extend([1, 1, 1, 1])

    out = g.node.add(name="classification", op="Identity", input=["conv_out"])
    out.attr["T"].type = types_pb2.DT_FLOAT

    # SignatureDef (serving_default) with a 7-channel input and a 3-class output.
    sig = meta.signature_def["serving_default"]
    in_info = sig.inputs["input_1"]
    in_info.name = "input_image:0"
    in_info.dtype = types_pb2.DT_FLOAT
    for d in (-1, 100, 221, 7):
        dim = in_info.tensor_shape.dim.add()
        dim.size = d

    out_info = sig.outputs["classification"]
    out_info.name = "classification:0"
    out_info.dtype = types_pb2.DT_FLOAT
    for d in (-1, 3):
        dim = out_info.tensor_shape.dim.add()
        dim.size = d

    # Write saved_model.pb
    out_dir = tmp / "synth_model"
    out_dir.mkdir(parents=True)
    (out_dir / "saved_model.pb").write_bytes(sm.SerializeToString())
    # A zeroed variables.index is a malformed tensor bundle — weights() should
    # raise when it tries to parse it (the magic/footer check fails).
    (out_dir / "variables").mkdir()
    (out_dir / "variables" / "variables.index").write_bytes(b"\x00" * 64)
    return out_dir


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        model = _build_synth_savedmodel(Path(tmpdir))
        reader = SavedModelReader(model)
        summary = reader.graph_summary()

        assert summary.signature_name == "serving_default", summary.signature_name
        assert len(summary.inputs) == 1, summary.inputs
        assert summary.inputs[0].name == "input_1"
        assert summary.inputs[0].dtype == "float32"
        assert summary.inputs[0].shape == [None, 100, 221, 7]
        assert len(summary.outputs) == 1
        assert summary.outputs[0].shape == [None, 3]

        ops = sorted({n.op for n in summary.nodes})
        assert ops == ["Conv2D", "Identity", "VarHandleOp"], ops

        assert "conv_kernel" in summary.variables, summary.variables

        # weights() is implemented now; on a zeroed/malformed variables.index
        # it must raise (the exact type depends on the TensorBundle reader —
        # e.g. RuntimeError on the magic check), but never succeed silently.
        try:
            reader.weights()
        except Exception as e:
            print(f"weights() correctly raises on a malformed index: {e!r}")
        else:
            print("FAIL: weights() should have raised on a zeroed variables.index")
            return 1

        print("OK")
        print(f"  signature: {summary.signature_name}")
        print(f"  inputs:    {summary.inputs}")
        print(f"  outputs:   {summary.outputs}")
        print(f"  nodes:     {len(summary.nodes)} ({ops})")
        print(f"  variables: {summary.variables}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
