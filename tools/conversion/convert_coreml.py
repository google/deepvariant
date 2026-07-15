"""Convert a DeepVariant WGS SavedModel to a Core ML .mlpackage — TF-free.

Pipeline:
  1. TensorBundle reader   (tensor_bundle_reader.py) reads weights.
  2. MIL program builder   (inception_v3_mil.py)     assembles the graph.
  3. coremltools.convert   wraps the MIL program into a .mlpackage.

Run inside venv-coreml (coremltools, no TF, no PyTorch):
    source tools/conversion/venv-coreml/bin/activate
    python convert_coreml.py \\
        --bundle models/wgs/variables/variables \\
        --output  models/wgs.mlpackage

Options:
    --compute-units   ALL (default, tries ANE first) | CPU_AND_GPU | CPU_ONLY
    --target          macOS14 (default) | macOS15
    --batch-min / --batch-max  (default 1 / 4096)
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time

try:
    import coremltools as ct
except ImportError as exc:
    raise ImportError(
        "coremltools not installed — activate venv-coreml first: "
        "source tools/conversion/venv-coreml/bin/activate"
    ) from exc

# Local imports — run from tools/conversion/ or add its path.
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "Generated"))

from tensor_bundle_reader import TensorBundle  # noqa: E402
import inception_v3_mil as iv3  # noqa: E402


_COMPUTE_UNITS = {
    "ALL": ct.ComputeUnit.ALL,
    "CPU_AND_GPU": ct.ComputeUnit.CPU_AND_GPU,
    "CPU_AND_NE": ct.ComputeUnit.CPU_AND_NE,
    "CPU_ONLY": ct.ComputeUnit.CPU_ONLY,
}

_TARGETS = {
    "macOS14": ct.target.macOS14,
    "macOS15": ct.target.macOS15,
}


def convert(
    bundle_prefix: str,
    output: str,
    compute_units: str = "ALL",
    target: str = "macOS14",
    batch_min: int = 1,
    batch_max: int = 4096,
) -> ct.models.MLModel:
    """Read weights, build MIL program, convert to Core ML .mlpackage."""
    print(f"loading bundle  {bundle_prefix}.{{index,data-*}}")
    t0 = time.time()
    bundle = TensorBundle(bundle_prefix)
    print(f"  {len(bundle.entries)} variables in {time.time()-t0:.1f}s")

    print("building MIL program …")
    t0 = time.time()
    prog = iv3.build_program(bundle, batch_min=batch_min, batch_max=batch_max)
    print(f"  done in {time.time()-t0:.1f}s")

    print(
        f"coremltools.convert  "
        f"(compute_units={compute_units}, target={target})"
    )
    t0 = time.time()
    dyn_shape = ct.Shape(
        shape=(
            ct.RangeDim(
                lower_bound=batch_min, upper_bound=batch_max, default=1,
            ),
            100, 221, 7,
        )
    )
    mlmodel = ct.convert(
        prog,
        inputs=[ct.TensorType(name="x", shape=dyn_shape)],
        compute_units=_COMPUTE_UNITS[compute_units],
        minimum_deployment_target=_TARGETS[target],
        convert_to="mlprogram",
    )
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s")

    mlmodel.author = "DeepVariant Apple Silicon Native Port"
    mlmodel.license = "BSD-3-Clause"
    mlmodel.short_description = (
        f"DeepVariant Inception-v3 (7-channel pileup, 3-class softmax) "
        f"built from WGS checkpoint via TF-free MIL path on "
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}"
    )

    if os.path.exists(output):
        print(f"removing existing  {output}")
        if os.path.isdir(output):
            shutil.rmtree(output)
        else:
            os.unlink(output)

    mlmodel.save(output)
    print(f"wrote  {output}")
    return mlmodel


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--bundle",
        required=True,
        help="TensorBundle prefix, e.g. models/wgs/variables/variables",
    )
    p.add_argument("--output", required=True, help="output .mlpackage path")
    p.add_argument(
        "--compute-units",
        default="ALL",
        choices=list(_COMPUTE_UNITS),
    )
    p.add_argument("--target", default="macOS14", choices=list(_TARGETS))
    p.add_argument("--batch-min", type=int, default=1)
    p.add_argument("--batch-max", type=int, default=4096)
    args = p.parse_args()
    convert(
        args.bundle,
        args.output,
        compute_units=args.compute_units,
        target=args.target,
        batch_min=args.batch_min,
        batch_max=args.batch_max,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
