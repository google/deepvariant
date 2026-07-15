"""Compare candidate softmax TFRecords to a reference, emit a diff report.

Usage:
    python parity_check.py \
        --reference ../reference/output/wgs/call_variants_chr20.tfrecord \
        --candidates ../../benchmarks/coreml_wgs.cv.tfrecord \
                     ../../benchmarks/metal_wgs.cv.tfrecord
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np


def _read_tfrecord_payloads(path: str):
    with open(path, "rb") as f:
        while True:
            length_bytes = f.read(8)
            if not length_bytes:
                return
            if len(length_bytes) != 8:
                raise RuntimeError(f"truncated tfrecord {path}")
            (length,) = struct.unpack("<Q", length_bytes)
            f.read(4)  # length crc — not verified in this dev tool
            payload = f.read(length)
            f.read(4)  # payload crc
            yield payload


def _decode_minimal(payload: bytes) -> tuple[int, np.ndarray]:
    """Decode the minimal {idx:u32, softmax:[3]f32} record written by bench.py."""
    if len(payload) != 4 + 12:
        raise RuntimeError(f"expected 16-byte payload, got {len(payload)}")
    idx = struct.unpack("<I", payload[:4])[0]
    sm = np.frombuffer(payload[4:], dtype=np.float32)
    return idx, sm


def load_softmax_minimal(path: str) -> np.ndarray:
    rows: dict[int, np.ndarray] = {}
    for payload in _read_tfrecord_payloads(path):
        idx, sm = _decode_minimal(payload)
        rows[idx] = sm
    if not rows:
        raise RuntimeError(f"no records in {path}")
    return np.stack([rows[i] for i in sorted(rows)], axis=0)


def load_softmax_dv(path: str) -> np.ndarray:
    """Decode upstream's CallVariantsOutput proto to extract the genotype probabilities.

    CallVariantsOutput.genotype_probabilities is a repeated double field.
    We use a minimal protobuf decoder rather than depending on the generated
    Python module — this script must work outside the conversion venv.
    """
    rows: list[np.ndarray] = []
    for payload in _read_tfrecord_payloads(path):
        # Walk the wire format looking for field 3 (genotype_probabilities, varint-tag 0x1A).
        # This is intentionally minimal — fragile if the proto schema changes.
        i = 0
        probs: list[float] = []
        while i < len(payload):
            tag = payload[i]
            i += 1
            field_no = tag >> 3
            wire = tag & 0x7
            if wire == 0:  # varint
                while payload[i] & 0x80:
                    i += 1
                i += 1
            elif wire == 1:  # 64-bit
                i += 8
            elif wire == 2:  # length-delimited
                ln = 0
                shift = 0
                while True:
                    b = payload[i]
                    i += 1
                    ln |= (b & 0x7F) << shift
                    if not (b & 0x80):
                        break
                    shift += 7
                seg = payload[i : i + ln]
                i += ln
                if field_no == 3:  # genotype_probabilities, packed doubles
                    probs.extend(np.frombuffer(seg, dtype=np.float64).tolist())
            elif wire == 5:  # 32-bit
                i += 4
            else:
                raise RuntimeError(f"unsupported wire type {wire}")
        if probs:
            rows.append(np.asarray(probs[:3], dtype=np.float32))
    if not rows:
        raise RuntimeError(f"no CallVariantsOutput records in {path}")
    return np.stack(rows, axis=0)


def compare(reference: np.ndarray, candidate: np.ndarray) -> dict:
    if reference.shape != candidate.shape:
        return {
            "ok": False,
            "reason": f"shape mismatch: ref {reference.shape} vs cand {candidate.shape}",
        }
    diff = np.abs(reference - candidate)
    max_abs = float(diff.max())
    mean_abs = float(diff.mean())
    ref_arg = reference.argmax(axis=1)
    cand_arg = candidate.argmax(axis=1)
    arg_disagree = int((ref_arg != cand_arg).sum())
    return {
        "ok": arg_disagree == 0 and max_abs <= 1e-3,
        "n": int(reference.shape[0]),
        "max_abs_softmax": max_abs,
        "mean_abs_softmax": mean_abs,
        "argmax_disagreements": arg_disagree,
        "argmax_disagreement_rate": arg_disagree / reference.shape[0],
    }


def autoload(path: str) -> np.ndarray:
    """Try the upstream CallVariantsOutput format first, fall back to minimal.

    The DV decoder is preferred. We only switch to the minimal decoder when
    the DV decode raises or yields zero rows, and we log which decoder won to
    stderr so a silent format switch never goes unnoticed.
    """
    try:
        rows = load_softmax_dv(path)
    except Exception as exc:
        print(
            f"parity_check: DV decode of {path} failed ({exc}); "
            f"falling back to minimal decoder",
            file=sys.stderr,
        )
        return load_softmax_minimal(path)
    if rows.shape[0] == 0:
        print(
            f"parity_check: DV decode of {path} returned zero rows; "
            f"falling back to minimal decoder",
            file=sys.stderr,
        )
        return load_softmax_minimal(path)
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True)
    p.add_argument("--candidates", nargs="+", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--max-abs-tol", type=float, default=1e-3)
    args = p.parse_args()

    print(f"reference: {args.reference}")
    ref = autoload(args.reference)
    print(f"  shape={ref.shape}")

    results: dict[str, dict] = {}
    overall_ok = True
    for c in args.candidates:
        print(f"candidate: {c}")
        cand = autoload(c)
        report = compare(ref, cand)
        results[c] = report
        if "reason" in report:
            # compare() bailed out early (e.g. shape mismatch) and the numeric
            # diff keys are absent — surface the reason instead of indexing them.
            report["ok"] = False
            print(f"  FAIL: {report['reason']}")
        else:
            report["ok"] = (
                report.get("argmax_disagreements", 1) == 0
                and report.get("max_abs_softmax", 1.0) <= args.max_abs_tol
            )
            line = (
                f"  n={report.get('n', '?')} max|Δ|={report.get('max_abs_softmax', 0):.2e} "
                f"mean|Δ|={report.get('mean_abs_softmax', 0):.2e} "
                f"argmax_disagree={report.get('argmax_disagreements', 0)}/{report.get('n', 0)} "
                f"({100 * report.get('argmax_disagreement_rate', 0):.3f}%) "
                f"{'OK' if report['ok'] else 'FAIL'}"
            )
            print(line)
        overall_ok &= report["ok"]

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump({"reference": args.reference, "results": results}, f, indent=2)
        print(f"wrote {args.output}")

    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
