"""SavedModel weights -> MLX-friendly safetensors. TF-free.

Reads the SavedModel via savedmodel_reader (pure protobuf, no TF), strips
the Keras name suffixes, and writes a safetensors bundle. The MLX
Inception-v3 architecture is rebuilt at bench time in
`bench.py --backend mlx`.

Run inside `venv-mlx` (MLX + safetensors, no TF).

Usage:
    python convert_mlx.py \\
        --saved-model models/wgs --output models/wgs.mlx.safetensors

Status: STUB. Depends on savedmodel_reader.py being implemented first.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

from savedmodel_reader import SavedModelReader


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--saved-model", required=True)
    p.add_argument("--output", required=True, help="output safetensors file")
    args = p.parse_args()

    reader = SavedModelReader(args.saved_model)
    try:
        weights = reader.weights()
    except NotImplementedError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    from safetensors.numpy import save_file  # noqa: WPS433

    metadata = {
        "captured_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "saved_model": os.path.abspath(args.saved_model),
        "n_tensors": str(len(weights)),
    }
    save_file(weights, args.output, metadata=metadata)
    sidecar = args.output + ".manifest.json"
    with open(sidecar, "w") as f:
        json.dump(
            {
                "metadata": metadata,
                "tensors": {
                    k: {"shape": list(v.shape), "dtype": str(v.dtype)}
                    for k, v in weights.items()
                },
            },
            f,
            indent=2,
        )
    print(f"wrote {args.output} ({len(weights)} tensors) and {sidecar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
