"""Dump per-layer reference outputs for a DeepVariant SavedModel.

Runs **inside google/deepvariant:1.10.0 Docker** (TF 2.16). Approach:

  1. Load the SavedModel and get the `serving_default` signature.
  2. Freeze the graph via `convert_variables_to_constants_v2` — this
     inlines the StatefulPartitionedCall function body into the outer
     graph and bakes all variables as Const nodes, so every Inception
     op becomes addressable by name at the top level (e.g.,
     `StatefulPartitionedCall/inceptionv3/mixed0/concat:0`).
  3. Run the frozen graph in a v1 `tf.Session` and fetch the input
     tensor + 19 named tap tensors + the final softmax in a single
     `session.run` call.
  4. Save every fetched tensor as `<tap>.npy` under `<out_dir>/`.
     For 4-D tensors the layout is converted from NHWC (TF native) to
     NCHW (Metal builder native) before saving.

Output (under `<out_dir>/`):
    _input.npy             fixed seed-0 input batch (NHWC)
    _savedmodel_softmax.npy   softmax via SavedModel sig (gold ref)
    _frozen_softmax.npy    softmax via frozen+v1 path (sanity check —
                           should match _savedmodel_softmax exactly)
    stem_s1a.npy ... gap.npy   per-tap NCHW FP32 dumps

The tap → tensor-name table is the **authoritative** mapping derived
by inspecting the frozen graph in the upstream Docker — no longer
relying on hand-coded `inception_v3_mil.py` indices.

Usage (inside Docker, see `dump_tf_per_layer.sh`):

    python3 dump_tf_per_layer.py /in /out
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf
import tensorflow.compat.v1 as tfv1
from tensorflow.python.framework.convert_to_constants import (
    convert_variables_to_constants_v2,
)


# ---------------------------------------------------------------------------
# Authoritative tap → frozen-graph tensor name mapping.
# Verified inside google/deepvariant:1.10.0 Docker by walking the inner
# StatefulPartitionedCall function body. See PORT_LOG.md.
# ---------------------------------------------------------------------------

_PFX = "StatefulPartitionedCall/inceptionv3"
TAPS: dict[str, str] = {
    "stem_s1a":  f"{_PFX}/activation/Relu:0",
    "stem_s2a":  f"{_PFX}/activation_1/Relu:0",
    "stem_s2b":  f"{_PFX}/activation_2/Relu:0",
    "stem_mp3a": f"{_PFX}/max_pooling2d/MaxPool:0",
    "stem_s3b":  f"{_PFX}/activation_3/Relu:0",
    "stem_s4a":  f"{_PFX}/activation_4/Relu:0",
    "stem_mp5a": f"{_PFX}/max_pooling2d_1/MaxPool:0",
    "5b":        f"{_PFX}/mixed0/concat:0",
    "5c":        f"{_PFX}/mixed1/concat:0",
    "5d":        f"{_PFX}/mixed2/concat:0",
    "6a":        f"{_PFX}/mixed3/concat:0",
    "6b":        f"{_PFX}/mixed4/concat:0",
    "6c":        f"{_PFX}/mixed5/concat:0",
    "6d":        f"{_PFX}/mixed6/concat:0",
    "6e":        f"{_PFX}/mixed7/concat:0",
    "7a":        f"{_PFX}/mixed8/concat:0",
    "7b":        f"{_PFX}/mixed9/concat:0",
    "7c":        f"{_PFX}/mixed10/concat:0",
    "gap":       f"{_PFX}/global_average_pooling2d/Mean:0",
}

# Input/output names of the frozen function.
INPUT_TENSOR = "input_1:0"
OUTPUT_TENSOR = "Identity:0"  # final softmax wrapped as Identity


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    if len(argv) not in (3, 4):
        print(
            f"usage: {argv[0]} <savedmodel_dir> <out_dir> [input.npy]\n"
            "  Without [input.npy], a fixed seed-0 random batch is used.\n"
            "  With [input.npy], it must be FP32 NHWC of shape (B, H, W, C)\n"
            "  matching the model's expected input geometry.",
            file=sys.stderr,
        )
        return 2

    model_dir = Path(argv[1])
    out_dir = Path(argv[2])
    in_npy = Path(argv[3]) if len(argv) == 4 else None
    out_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect input shape from upstream's example_info.json.
    info = model_dir / "model.example_info.json"
    if info.exists():
        sh = json.loads(info.read_text())["shape"]
        H, W, C = int(sh[0]), int(sh[1]), int(sh[2])
    else:
        H, W, C = 100, 221, 7

    # Either load real pileup batch (preferred for drift profiling) or
    # fall back to fixed-seed random.
    if in_npy is not None:
        x = np.load(in_npy).astype(np.float32)
        if x.ndim != 4 or x.shape[1:] != (H, W, C):
            print(
                f"input.npy shape mismatch: got {x.shape}, "
                f"expected (B, {H}, {W}, {C})",
                file=sys.stderr,
            )
            return 2
        print(f"loaded input from {in_npy}: shape={x.shape}")
    else:
        rng = np.random.default_rng(0)
        x = rng.uniform(0.0, 255.0, (1, H, W, C)).astype(np.float32)
        print(f"input shape: {x.shape} (seed-0 random)")
    np.save(out_dir / "_input.npy", x)

    # 1) Canonical SavedModel signature forward (gold reference).
    print(f"loading SavedModel: {model_dir}")
    sm = tf.saved_model.load(str(model_dir))
    fn = sm.signatures["serving_default"]
    sm_out = fn(input_1=tf.constant(x))
    sm_softmax = list(sm_out.values())[0].numpy()
    np.save(out_dir / "_savedmodel_softmax.npy", sm_softmax)
    print(f"  SavedModel softmax: {sm_softmax[0]}")

    # 2) Freeze the graph (inlines functions + bakes variables to Const).
    print("freezing graph (convert_variables_to_constants_v2) ...")
    frozen = convert_variables_to_constants_v2(fn)
    gd = frozen.graph.as_graph_def()
    print(
        f"  frozen graph: {len(gd.node)} top-level ops, "
        f"{len(gd.library.function)} fns"
    )

    # 3) Re-import the frozen GraphDef into a v1 Graph, run with fetches.
    print("running v1 Session with intermediate fetches ...")
    g = tfv1.Graph()
    with g.as_default():
        tfv1.import_graph_def(gd, name="")

    fetches = list(TAPS.values()) + [OUTPUT_TENSOR]
    with tfv1.Session(graph=g) as sess:
        results = sess.run(fetches, feed_dict={INPUT_TENSOR: x})

    *tap_arrays, frozen_softmax = results

    # 4) Save softmax sanity check.
    np.save(out_dir / "_frozen_softmax.npy", frozen_softmax)
    diff = float(np.abs(sm_softmax - frozen_softmax).max())
    print(f"  Frozen softmax: {frozen_softmax[0]}")
    print(f"  max-abs diff (frozen vs SavedModel sig): {diff:.6e}")
    if diff > 1e-4:
        print(
            "WARN: frozen graph diverged from SavedModel signature — "
            "freezing introduced a numerical regression."
        )

    # 5) Save each tap (NHWC for 4-D — matches Metal builder's native
    # layout — and identity for 2-D).
    print("saving taps ...")
    for (tap_name, _tensor_name), arr in zip(TAPS.items(), tap_arrays):
        arr_save = arr.astype(np.float32)
        np.save(out_dir / f"{tap_name}.npy", arr_save)
        print(
            f"  {tap_name:<11} shape={arr_save.shape} "
            f"({arr_save.nbytes // 1024} KB)"
        )

    print(f"\nsaved {len(TAPS)} taps to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
