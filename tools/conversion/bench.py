"""Bench inference latency, throughput, and GPU/ANE residency — TF-free.

Reads input examples from a TFRecord (raw protobuf, no TF runtime), feeds
them into the chosen backend (Core ML or MLX), captures softmax outputs and
powermetrics GPU/ANE residency, writes one JSON metrics file and one
TFRecord of softmax records (for parity_check.py).

Usage:
    python bench.py --backend coreml \\
        --model models/wgs.mlpackage \\
        --examples ../reference/cache/wgs_chr20_1000.tfrecord \\
        --output ../../benchmarks/coreml_wgs.json \\
        --output-cv ../../benchmarks/coreml_wgs.cv.tfrecord
"""

from __future__ import annotations

import argparse
import json
import signal
import struct
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator

import numpy as np


# ---------------------------------------------------------------------------
# TFRecord I/O — raw, no TF
# ---------------------------------------------------------------------------

_CRC_MASK = 0xA282EAD8


def _crc32c(data: bytes) -> int:
    import google_crc32c

    c = google_crc32c.Checksum()
    c.update(data)
    return int.from_bytes(c.digest(), "big")


def _masked_crc32(data: bytes) -> int:
    crc = _crc32c(data)
    return ((crc >> 15) | ((crc << 17) & 0xFFFFFFFF)) + _CRC_MASK & 0xFFFFFFFF


def write_tfrecord(path: str, payloads: Iterator[bytes]) -> int:
    n = 0
    with open(path, "wb") as f:
        for payload in payloads:
            length = struct.pack("<Q", len(payload))
            f.write(length)
            f.write(struct.pack("<I", _masked_crc32(length)))
            f.write(payload)
            f.write(struct.pack("<I", _masked_crc32(payload)))
            n += 1
    return n


def read_tfrecord(path: str) -> Iterator[bytes]:
    """Yield each record's payload bytes. CRCs are not verified (dev tool)."""
    with open(path, "rb") as f:
        while True:
            ln_b = f.read(8)
            if not ln_b:
                return
            if len(ln_b) != 8:
                raise RuntimeError(f"truncated tfrecord {path}")
            (length,) = struct.unpack("<Q", ln_b)
            f.read(4)
            payload = f.read(length)
            f.read(4)
            yield payload


# ---------------------------------------------------------------------------
# Minimal protobuf decoder for tf.train.Example
#
# Once we vendor TF .proto files and run protoc (see Protos/README.md), this
# block is replaced by `from tensorflow_pb2 import Example`. Until then this
# hand-rolled walker handles the small subset of fields we actually consume.
# ---------------------------------------------------------------------------

def _read_varint(buf: bytes, i: int) -> tuple[int, int]:
    val = 0
    shift = 0
    while True:
        b = buf[i]
        i += 1
        val |= (b & 0x7F) << shift
        if not (b & 0x80):
            return val, i
        shift += 7


def parse_tf_example(payload: bytes, h: int, w: int, c: int) -> np.ndarray:
    """Parse a tf.train.Example proto and return its image as (H,W,C) float32."""
    target_keys = {"image/encoded", "image"}
    i = 0
    while i < len(payload):
        tag = payload[i]
        i += 1
        field, wire = tag >> 3, tag & 7
        if wire == 2:
            ln, i = _read_varint(payload, i)
            seg = payload[i:i + ln]
            i += ln
            if field == 1:  # Example.features
                arr = _scan_features(seg, target_keys, h, w, c)
                if arr is not None:
                    return arr
        elif wire == 0:
            _, i = _read_varint(payload, i)
        elif wire == 1:
            i += 8
        elif wire == 5:
            i += 4
        else:
            raise RuntimeError(f"unsupported wire type {wire}")
    raise RuntimeError("no image feature in example")


def _scan_features(
    seg: bytes,
    keys: set[str],
    h: int,
    w: int,
    c: int,
) -> np.ndarray | None:
    """Walk Features { map<string,Feature> feature = 1; } looking for `keys`."""
    i = 0
    while i < len(seg):
        tag = seg[i]
        i += 1
        if tag & 7 != 2:
            raise RuntimeError("Features map entries must be length-delimited")
        ln, i = _read_varint(seg, i)
        entry = seg[i:i + ln]
        i += ln
        key, value = None, None
        j = 0
        while j < len(entry):
            t = entry[j]
            j += 1
            f, wire = t >> 3, t & 7
            if wire != 2:
                raise RuntimeError("MapEntry fields are length-delimited")
            l2, j = _read_varint(entry, j)
            buf = entry[j:j + l2]
            j += l2
            if f == 1:
                key = buf.decode("utf-8")
            elif f == 2:
                value = buf
        if key in keys and value is not None:
            return _decode_feature(
                value, h, w, c, prefer_bytes=(key == "image/encoded"),
            )
    return None


def _decode_feature(
    buf: bytes,
    h: int,
    w: int,
    c: int,
    prefer_bytes: bool,
) -> np.ndarray:
    """Feature is a oneof of {bytes_list=1, float_list=2, int64_list=3}."""
    i = 0
    while i < len(buf):
        tag = buf[i]
        i += 1
        field, wire = tag >> 3, tag & 7
        if wire != 2:
            raise RuntimeError("Feature list fields are length-delimited")
        ln, i = _read_varint(buf, i)
        seg = buf[i:i + ln]
        i += ln
        if field == 1 and prefer_bytes:  # BytesList
            j = 1  # skip BytesList.value tag
            l2, j = _read_varint(seg, j)
            data = seg[j:j + l2]
            pixels = np.frombuffer(data, dtype=np.uint8)
            if pixels.size != h * w * c:
                raise RuntimeError(
                    f"image/encoded BytesList has {pixels.size} uint8 elements, "
                    f"expected {h * w * c} for shape ({h}, {w}, {c}); the feature "
                    "may be a compressed encoding (e.g. PNG) rather than raw "
                    "uint8 pixels"
                )
            return pixels.astype(np.float32).reshape(h, w, c)
        if field == 2 and not prefer_bytes:  # FloatList (packed)
            floats = np.frombuffer(seg, dtype=np.float32)
            if floats.size != h * w * c:
                raise RuntimeError(
                    f"FloatList has {floats.size} float32 elements, expected "
                    f"{h * w * c} for shape ({h}, {w}, {c})"
                )
            return floats.reshape(h, w, c)
    raise RuntimeError("Feature has no recognized list")


def read_examples(path: str, input_shape: tuple[int, int, int]) -> np.ndarray:
    h, w, c = input_shape
    arrs = [parse_tf_example(p, h, w, c) for p in read_tfrecord(path)]
    if not arrs:
        raise RuntimeError(f"no examples in {path}")
    return np.stack(arrs, axis=0)


# ---------------------------------------------------------------------------
# powermetrics side-thread — captures GPU and ANE power residency.
# ---------------------------------------------------------------------------

@dataclass
class GpuStats:
    samples: int
    gpu_power_mean_mw: float
    gpu_power_max_mw: float
    ane_power_mean_mw: float
    ane_power_max_mw: float


class PowerSampler:
    """Run powermetrics in background. Requires sudo (-n)."""

    def __init__(self, interval_ms: int = 500) -> None:
        self.interval_ms = interval_ms
        self.proc: subprocess.Popen | None = None
        self.samples: list[tuple[float, float]] = []
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def __enter__(self) -> "PowerSampler":
        cmd = [
            "sudo", "-n", "powermetrics",
            "--samplers", "gpu_power,ane_power",
            "-i", str(self.interval_ms),
            "-f", "text",
        ]
        try:
            self.proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
            )
        except FileNotFoundError:
            print("warning: powermetrics not found", file=sys.stderr)
            self.proc = None
            return self
        self._thread = threading.Thread(target=self._reader, daemon=True)
        self._thread.start()
        return self

    def _reader(self) -> None:
        assert self.proc is not None
        block: list[str] = []
        for line in self.proc.stdout:  # type: ignore[union-attr]
            if self._stop.is_set():
                break
            block.append(line)
            if line.startswith("**") and len(block) > 5:
                gpu = ane = 0.0
                for ln in block:
                    if ln.startswith("GPU Power:"):
                        try:
                            gpu = float(ln.split(":", 1)[1].strip().split()[0])
                        except Exception:
                            pass
                    elif ln.startswith("ANE Power:"):
                        try:
                            ane = float(ln.split(":", 1)[1].strip().split()[0])
                        except Exception:
                            pass
                if gpu or ane:
                    self.samples.append((gpu, ane))
                block = []

    def __exit__(self, *exc) -> None:
        self._stop.set()
        if self.proc is not None:
            try:
                self.proc.send_signal(signal.SIGINT)
                self.proc.wait(timeout=2)
            except Exception:
                self.proc.kill()
        if self._thread is not None:
            self._thread.join(timeout=2)

    def stats(self) -> GpuStats:
        if not self.samples:
            return GpuStats(0, 0.0, 0.0, 0.0, 0.0)
        g = np.array([s[0] for s in self.samples])
        a = np.array([s[1] for s in self.samples])
        return GpuStats(
            samples=len(self.samples),
            gpu_power_mean_mw=float(g.mean()),
            gpu_power_max_mw=float(g.max()),
            ane_power_mean_mw=float(a.mean()),
            ane_power_max_mw=float(a.max()),
        )


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

def bench_coreml(
    model_path: str, x: np.ndarray, batch: int,
) -> tuple[np.ndarray, float]:
    import coremltools as ct

    m = ct.models.MLModel(model_path, compute_units=ct.ComputeUnit.ALL)
    in_name = m.get_spec().description.input[0].name
    out_chunks: list[np.ndarray] = []
    t0 = time.time()
    for i in range(0, x.shape[0], batch):
        chunk = x[i:i + batch]
        result = m.predict({in_name: chunk})
        out_chunks.append(
            np.asarray(next(iter(result.values())), dtype=np.float32),
        )
    return np.concatenate(out_chunks, axis=0), time.time() - t0


def bench_mlx(
    model_path: str, x: np.ndarray, batch: int,
) -> tuple[np.ndarray, float]:
    raise NotImplementedError(
        "mlx bench: needs hand-built MLX Inception-v3 module loaded from "
        "the safetensors weights produced by convert_mlx.py. "
        "Phase 0 sub-step pending.",
    )


def write_softmax_records(path: str, softmax: np.ndarray) -> int:
    def gen() -> Iterator[bytes]:
        for i, row in enumerate(softmax):
            yield struct.pack("<I", i) + row.astype(np.float32).tobytes()

    return write_tfrecord(path, gen())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", required=True, choices=["coreml", "mlx"])
    p.add_argument("--model", required=True)
    p.add_argument("--examples", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--output-cv", required=True)
    p.add_argument("--batch", type=int, default=128)
    p.add_argument("--input-shape", default="100,221,7")
    p.add_argument("--warmup-batches", type=int, default=2)
    p.add_argument("--no-powermetrics", action="store_true")
    args = p.parse_args()

    h, w, c = (int(s) for s in args.input_shape.split(","))
    print(f"loading examples from {args.examples}")
    x = read_examples(args.examples, (h, w, c))
    print(f"  shape={x.shape}, dtype={x.dtype}")

    print(f"warmup {args.warmup_batches} batches @ batch={args.batch}")
    warm_n = min(args.warmup_batches * args.batch, x.shape[0])
    if args.backend == "coreml":
        bench_coreml(args.model, x[:warm_n], args.batch)
    elif args.backend == "mlx":
        bench_mlx(args.model, x[:warm_n], args.batch)

    print("benching...")
    sampler = PowerSampler() if not args.no_powermetrics else None
    if sampler is not None:
        sampler.__enter__()
    try:
        if args.backend == "coreml":
            softmax, elapsed = bench_coreml(args.model, x, args.batch)
        elif args.backend == "mlx":
            softmax, elapsed = bench_mlx(args.model, x, args.batch)
    finally:
        if sampler is not None:
            sampler.__exit__(None, None, None)

    n = x.shape[0]
    print(f"  {n} examples in {elapsed:.3f}s = {n / elapsed:.1f} ex/s")

    Path(args.output_cv).parent.mkdir(parents=True, exist_ok=True)
    n_written = write_softmax_records(args.output_cv, softmax)
    print(f"  wrote {n_written} softmax records to {args.output_cv}")

    metrics = {
        "backend": args.backend,
        "model": args.model,
        "examples": args.examples,
        "n_examples": n,
        "batch_size": args.batch,
        "elapsed_seconds": elapsed,
        "examples_per_second": n / elapsed,
        "input_shape": [h, w, c],
        "softmax_shape": list(softmax.shape),
        "softmax_dtype": str(softmax.dtype),
    }
    if sampler is not None:
        metrics["gpu_stats"] = asdict(sampler.stats())

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"wrote metrics to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
