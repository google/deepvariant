#!/usr/bin/env python3
"""Pack a TF SavedModel TensorBundle into a single .dvw weight file.

`.dvw` is the runtime weight format consumed by the Phase 5.5 Metal/BNNS
inference path. Single file, deterministic byte layout, no TF runtime
needed at load time. Only DT_FLOAT (dtype=1) tensors are kept (model
weights + biases + BN params); the `_CHECKPOINTABLE_OBJECT_GRAPH`
serialized graph blob (dtype=7 STRING) is dropped.

Layout (all integers little-endian):

    +--------------------------------------------------------------+
    | magic[4]    = 'DVW1'                                         |
    | version[4]  = 1                                              |
    | n_tensors[4]                                                 |
    +--------------------------------------------------------------+
    | for each tensor (sorted by name for determinism):            |
    |   name_len[4]                                                |
    |   name[name_len] (utf-8)                                     |
    |   dtype[1] = 1 (DT_FLOAT only for now)                       |
    |   ndim[1]                                                    |
    |   shape[ndim*4] (uint32 le, in source/HWIO order)            |
    |   offset[8] (into payload, after the table)                  |
    |   n_bytes[8]                                                 |
    +--------------------------------------------------------------+
    | payload: concatenated raw FP32 LE bytes, in tensor order     |
    +--------------------------------------------------------------+

Usage:
    extract_weights.py <savedmodel_dir> <output.dvw>

Reads <savedmodel_dir>/variables/variables.{index, data-NNNNN-of-MMMMM}
via the existing `tensor_bundle_reader` (TF-free path through the
vendored protos in Generated/).
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

# tensor_bundle_reader prepends Generated/ to sys.path itself, which is what
# pulls in the TF-free vendored proto bindings.
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE / "Generated"))

import numpy as np

from tensor_bundle_reader import TensorBundle  # noqa: E402

DVW_MAGIC = b"DVW1"
DVW_VERSION = 1
DT_FLOAT = 1  # tensorflow DataType enum


def pack_dvw(savedmodel_dir: Path, out_path: Path) -> None:
    bundle_prefix = savedmodel_dir / "variables" / "variables"
    bundle = TensorBundle(bundle_prefix)

    # Keep only float weights, sorted by name for byte-stable output.
    keep = sorted(
        n for n, e in bundle.entries.items() if e.dtype == DT_FLOAT
    )
    if not keep:
        raise RuntimeError(f"no DT_FLOAT tensors found in {bundle_prefix}")

    # First pass: assemble payload + per-tensor offset/size, gather bytes
    # in a list to avoid materialising the 80+ MB blob twice.
    payload_chunks: list[bytes] = []
    cursor = 0
    metadata: list[tuple[str, np.ndarray, int, int]] = []
    for name in keep:
        arr = bundle.read_tensor(name).astype(np.float32, copy=False)
        # Force little-endian. M-series is little-endian natively but be
        # explicit so the file is portable to any consumer.
        if arr.dtype.byteorder == ">":
            arr = arr.astype("<f4")
        raw = arr.tobytes(order="C")
        metadata.append((name, arr, cursor, len(raw)))
        payload_chunks.append(raw)
        cursor += len(raw)

    # Second pass: write header + per-tensor table + payload.
    with open(out_path, "wb") as f:
        f.write(DVW_MAGIC)
        f.write(struct.pack("<I", DVW_VERSION))
        f.write(struct.pack("<I", len(metadata)))
        for name, arr, offset, n_bytes in metadata:
            name_bytes = name.encode("utf-8")
            f.write(struct.pack("<I", len(name_bytes)))
            f.write(name_bytes)
            f.write(struct.pack("<B", DT_FLOAT))
            f.write(struct.pack("<B", arr.ndim))
            for d in arr.shape:
                f.write(struct.pack("<I", int(d)))
            f.write(struct.pack("<Q", offset))
            f.write(struct.pack("<Q", n_bytes))
        for chunk in payload_chunks:
            f.write(chunk)

    n = len(metadata)
    total = cursor
    print(f"wrote {out_path} ({n} tensors, {total} payload bytes)")


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: extract_weights.py <savedmodel_dir> <out.dvw>",
              file=sys.stderr)
        return 2
    pack_dvw(Path(argv[1]), Path(argv[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
