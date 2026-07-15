"""TF-free SavedModel reader.

Two halves:

1. **Graph parsing** (this file) — implemented. Uses the protobuf bindings
   generated from `Protos/tensorflow/` under `Generated/` to read
   `saved_model.pb` and extract:
       - SignatureDef (input/output names + dtypes + shapes)
       - GraphDef.node[] (op type, attrs, inputs)
       - List of variable names referenced by the graph

2. **Weight extraction** (TensorBundle reader, follow-up commit) — STUB.
   Needs a small SSTable parser to read the `variables/variables.{index,
   data-00000-of-00001}` pair. Format is documented at
   `tensorflow/core/util/tensor_bundle/tensor_bundle.h`. The index is an
   LSM-style sstable mapping `variable_name -> BundleEntryProto`; each
   entry has `(offset, size, dtype, shape, slices)` into the data shard.

Importing this module **requires** that `Generated/` is populated:
    cd tools/conversion
    rm -rf Generated && mkdir Generated
    protoc --python_out=Generated -I=Protos \\
        $(cd Protos && find tensorflow -name '*.proto')
    touch Generated/__init__.py

(setup_venvs.sh runs this automatically.)
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Make the generated bindings importable regardless of where this script
# is invoked from. We deliberately don't pip-install them; they're regen
# every venv setup.
_GEN = Path(__file__).resolve().parent / "Generated"
if str(_GEN) not in sys.path:
    sys.path.insert(0, str(_GEN))


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

@dataclass
class TensorSpec:
    name: str
    dtype: str
    shape: list[int | None]


@dataclass
class NodeInfo:
    name: str
    op: str
    inputs: list[str]
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass
class GraphSummary:
    """A high-level view of one MetaGraphDef.

    `signature` covers the user-facing serving inputs/outputs
    (`serving_default`). `nodes` is the full op list. `variables` is the
    names of variables the graph references — what TensorBundle keys
    against.
    """
    signature_name: str
    inputs: list[TensorSpec]
    outputs: list[TensorSpec]
    nodes: list[NodeInfo]
    variables: list[str]


# ---------------------------------------------------------------------------
# SavedModel reader
# ---------------------------------------------------------------------------

class SavedModelReader:
    def __init__(self, directory: str | Path) -> None:
        self.directory = Path(directory)
        pb = self.directory / "saved_model.pb"
        if not pb.exists():
            raise FileNotFoundError(f"{pb} not found")

        # Lazy import — useful error if protoc bindings aren't generated.
        try:
            from tensorflow.core.protobuf import (  # type: ignore
                saved_model_pb2,
            )
        except ImportError as e:
            raise ImportError(
                "TF protobuf bindings not generated. Run setup_venvs.sh "
                "or the protoc invocation in Protos/SOURCES.md."
            ) from e

        sm = saved_model_pb2.SavedModel()
        with open(pb, "rb") as f:
            sm.ParseFromString(f.read())
        if not sm.meta_graphs:
            raise RuntimeError(f"{pb} has no MetaGraphDefs")
        self._sm = sm
        # Most DV SavedModels expose a single MetaGraphDef.
        self._meta = sm.meta_graphs[0]

    # ------------------------------------------------------------------ graph

    def graph_summary(
        self, signature: str = "serving_default",
    ) -> GraphSummary:
        sig_map = self._meta.signature_def
        if signature not in sig_map:
            available = list(sig_map.keys())
            raise KeyError(
                f"signature {signature!r} not found; available: {available}"
            )
        sig = sig_map[signature]

        inputs = [
            _tensor_info_to_spec(name, ti)
            for name, ti in sig.inputs.items()
        ]
        outputs = [
            _tensor_info_to_spec(name, ti)
            for name, ti in sig.outputs.items()
        ]

        graph_def = self._meta.graph_def
        nodes: list[NodeInfo] = []
        for n in graph_def.node:
            nodes.append(NodeInfo(
                name=n.name,
                op=n.op,
                inputs=list(n.input),
                attrs=_attrs_to_dict(n.attr),
            ))

        variables = sorted(_collect_variable_names(self._meta))

        return GraphSummary(
            signature_name=signature,
            inputs=inputs,
            outputs=outputs,
            nodes=nodes,
            variables=variables,
        )

    # ---------------------------------------------------------------- weights

    def weights(self) -> dict[str, Any]:
        """Read variable tensors from variables/variables.{index, data-*}.

        Returns a dict mapping variable name (e.g. ``conv2d_1/kernel``) to
        a numpy.ndarray with the original dtype + shape. No TF runtime —
        uses the in-tree TensorBundle reader.
        """
        idx = self.directory / "variables" / "variables.index"
        if not idx.exists():
            raise FileNotFoundError(f"{idx} not found")

        from tensor_bundle_reader import TensorBundle

        bundle = TensorBundle(self.directory / "variables" / "variables")
        return {name: bundle.read_tensor(name) for name in bundle.names()}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# DataType enum mapping (TF -> human-friendly).
_DTYPE_NAMES = {
    1: "float32", 2: "float64", 3: "int32", 4: "uint8", 5: "int16",
    6: "int8", 7: "string", 8: "complex64", 9: "int64", 10: "bool",
    14: "bfloat16", 15: "uint16", 16: "complex128", 17: "half",
    19: "uint32", 22: "uint64",
}


def _tensor_info_to_spec(name: str, ti) -> TensorSpec:
    """Convert a TensorInfo proto to our TensorSpec."""
    shape = []
    if ti.HasField("tensor_shape"):
        for dim in ti.tensor_shape.dim:
            shape.append(dim.size if dim.size >= 0 else None)
    return TensorSpec(
        name=name,
        dtype=_DTYPE_NAMES.get(ti.dtype, f"DT_{ti.dtype}"),
        shape=shape,
    )


def _attrs_to_dict(attrs) -> dict[str, Any]:
    """Convert a NodeDef.attr map to a plain dict.

    We only extract a few useful kinds (int, float, string, bool, type,
    shape, list of ints). Tensors are returned as a placeholder dict —
    don't materialise their bytes here, the bundle reader handles those.
    """
    out: dict[str, Any] = {}
    for key, value in attrs.items():
        kind = value.WhichOneof("value")
        if kind == "i":
            out[key] = int(value.i)
        elif kind == "f":
            out[key] = float(value.f)
        elif kind == "b":
            out[key] = bool(value.b)
        elif kind == "s":
            out[key] = value.s.decode("utf-8", errors="replace")
        elif kind == "type":
            out[key] = _DTYPE_NAMES.get(value.type, f"DT_{value.type}")
        elif kind == "shape":
            out[key] = [
                d.size if d.size >= 0 else None for d in value.shape.dim
            ]
        elif kind == "list":
            lst = value.list
            if lst.i:
                out[key] = list(lst.i)
            elif lst.f:
                out[key] = list(lst.f)
            elif lst.s:
                out[key] = [
                    s.decode("utf-8", errors="replace") for s in lst.s
                ]
            elif lst.type:
                out[key] = [
                    _DTYPE_NAMES.get(t, f"DT_{t}") for t in lst.type
                ]
            else:
                out[key] = "<list>"
        elif kind == "tensor":
            t = value.tensor
            shape = [
                d.size if d.size >= 0 else None for d in t.tensor_shape.dim
            ]
            out[key] = {
                "kind": "tensor",
                "dtype": _DTYPE_NAMES.get(t.dtype, f"DT_{t.dtype}"),
                "shape": shape,
            }
        else:
            out[key] = f"<{kind}>"
    return out


def _collect_variable_names(meta) -> set[str]:
    """Find variable names this MetaGraphDef references.

    For TF 2.x SavedModels: each variable is a `VarHandleOp` node whose
    `shared_name` attr is the bundle key. We collect those, plus the
    classic `variable_def` proto entries used by older bundles.
    """
    names: set[str] = set()
    for n in meta.graph_def.node:
        if n.op in ("VarHandleOp", "VariableV2", "VarHandleOp_v2"):
            attrs = _attrs_to_dict(n.attr)
            names.add(attrs.get("shared_name") or n.name)

    # Fallback: collection_def["variables"] holds serialized VariableDef.
    if "variables" in meta.collection_def:
        cd = meta.collection_def["variables"]
        for raw in cd.bytes_list.value:
            try:
                from tensorflow.core.framework import (  # type: ignore
                    variable_pb2,
                )
                vd = variable_pb2.VariableDef()
                vd.ParseFromString(raw)
                if vd.variable_name:
                    names.add(vd.variable_name)
            except Exception:
                pass
    return names
