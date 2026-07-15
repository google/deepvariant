"""Compare two CallVariantsOutput TFRecord files.

Reports:
  - Records counts on each side
  - Argmax agreement (% of variants where the called genotype matches)
  - Softmax max-abs diff (overall and per-variant 95th percentile)
  - List of disagreements (top by max-abs diff)

Usage:
    python parity_check.py REFERENCE_CVO_TFRECORD CANDIDATE_CVO_TFRECORD
"""

from __future__ import annotations

import struct
import sys
from collections import Counter


def _read_records(path: str):
    """Yield raw payload bytes from an uncompressed TFRecord."""
    with open(path, "rb") as f:
        while True:
            head = f.read(12)  # uint64 length + uint32 crc(length)
            if not head:
                return
            if len(head) != 12:
                raise RuntimeError(f"truncated tfrecord: {path}")
            (length,) = struct.unpack("<Q", head[:8])
            payload = f.read(length)
            f.read(4)  # payload crc — not validated here
            if len(payload) != length:
                raise RuntimeError(f"short payload in {path}")
            yield payload


def _read_varint(buf: bytes, i: int) -> tuple[int, int]:
    res = 0
    shift = 0
    while True:
        b = buf[i]
        i += 1
        res |= (b & 0x7F) << shift
        if not (b & 0x80):
            return res, i
        shift += 7


def _parse_cvo(payload: bytes) -> tuple[str, list[float]]:
    """Extract (variant_key, genotype_probs) from a CallVariantsOutput proto.

    variant_key = '<reference_name>:<start>:<ref>>:<alt>' for sorting/lookup.
    Probabilities are repeated double in field 3 (packed wire type 2 or
    one fixed64 per entry).
    """
    n = len(payload)
    i = 0
    variant_ref = ""
    variant_start = 0
    variant_ref_bases = ""
    variant_alts: list[str] = []
    probs: list[float] = []

    while i < n:
        tag, i = _read_varint(payload, i)
        field = tag >> 3
        wire = tag & 7
        if wire == 2:
            seg_len, i = _read_varint(payload, i)
            sub = payload[i:i + seg_len]
            i += seg_len
            if field == 1:  # variant message
                ref, start, refbases, alts = _parse_variant(sub)
                variant_ref, variant_start = ref, start
                variant_ref_bases, variant_alts = refbases, alts
            elif field == 3:  # packed double
                count = seg_len // 8
                probs = list(struct.unpack(f"<{count}d", sub))
            # Skip alt_allele_indices (field 2) and other fields.
        elif wire == 1:  # fixed64 → unpacked double
            d = struct.unpack("<d", payload[i:i + 8])[0]
            i += 8
            if field == 3:
                probs.append(d)
        elif wire == 0:
            _, i = _read_varint(payload, i)
        else:
            break

    alt_str = ",".join(variant_alts)
    key = f"{variant_ref}:{variant_start}:{variant_ref_bases}>{alt_str}"
    return key, probs


def _parse_variant(buf: bytes) -> tuple[str, int, str, list[str]]:
    """Extract reference_name, start, reference_bases, alternate_bases from
    a nucleus.genomics.v1.Variant payload (variants.proto field tags)."""
    n = len(buf)
    i = 0
    ref_name = ""
    start = 0
    ref_bases = ""
    alt_bases: list[str] = []
    while i < n:
        tag, i = _read_varint(buf, i)
        field = tag >> 3
        wire = tag & 7
        if wire == 2:
            seg_len, i = _read_varint(buf, i)
            sub = buf[i:i + seg_len].decode("utf-8", errors="replace")
            i += seg_len
            if field == 14:    # reference_name
                ref_name = sub
            elif field == 6:   # reference_bases
                ref_bases = sub
            elif field == 7:   # alternate_bases (repeated)
                alt_bases.append(sub)
        elif wire == 0:
            v, i = _read_varint(buf, i)
            if field == 16:    # start (int64)
                start = v
        elif wire == 1:        # fixed64 — skip
            i += 8
        elif wire == 5:        # fixed32 — skip
            i += 4
        else:
            break
    return ref_name, start, ref_bases, alt_bases


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2

    ref_path, cand_path = argv[1], argv[2]
    ref_records = {}
    for payload in _read_records(ref_path):
        key, probs = _parse_cvo(payload)
        ref_records[key] = probs

    cand_records = {}
    for payload in _read_records(cand_path):
        key, probs = _parse_cvo(payload)
        cand_records[key] = probs

    print(f"reference : {len(ref_records)} records  ({ref_path})")
    print(f"candidate : {len(cand_records)} records  ({cand_path})")

    common = sorted(set(ref_records.keys()) & set(cand_records.keys()))
    only_ref = sorted(set(ref_records.keys()) - set(cand_records.keys()))
    only_cand = sorted(set(cand_records.keys()) - set(ref_records.keys()))
    print(f"common    : {len(common)}")
    if only_ref:
        print(f"only in ref : {len(only_ref)} (e.g. {only_ref[:3]})")
    if only_cand:
        print(f"only in cand: {len(only_cand)} (e.g. {only_cand[:3]})")

    if not common:
        print("no overlapping variants — nothing to compare")
        return 1

    argmax_disagree = 0
    max_abs = 0.0
    abs_diffs: list[tuple[float, str, list[float], list[float]]] = []
    arg_dist = Counter()
    for k in common:
        rp = ref_records[k]
        cp = cand_records[k]
        if not rp or not cp or len(rp) != len(cp):
            argmax_disagree += 1
            continue
        ra = max(range(len(rp)), key=lambda j: rp[j])
        ca = max(range(len(cp)), key=lambda j: cp[j])
        if ra != ca:
            argmax_disagree += 1
        arg_dist[(ra, ca)] += 1
        d = max(abs(a - b) for a, b in zip(rp, cp))
        max_abs = max(max_abs, d)
        abs_diffs.append((d, k, rp, cp))

    n_common = len(common)
    pct = 100.0 * (n_common - argmax_disagree) / n_common
    print()
    n_match = n_common - argmax_disagree
    print(f"argmax agreement : {n_match}/{n_common} = {pct:.3f}%")
    print(f"softmax max-abs  : {max_abs:.6f}")

    abs_diffs.sort(reverse=True)
    if abs_diffs:
        print("\ntop 5 worst per-variant max-abs diffs:")
        for d, k, rp, cp in abs_diffs[:5]:
            print(f"  {d:.4f}  {k}")
            print(f"    ref :  {tuple(round(p, 4) for p in rp)}")
            print(f"    cand:  {tuple(round(p, 4) for p in cp)}")

    print("\n(ref_argmax, cand_argmax) -> count")
    for (ra, ca), c in sorted(arg_dist.items()):
        marker = "" if ra == ca else "  ← MISMATCH"
        print(f"  ({ra}, {ca})  {c}{marker}")

    if pct >= 100.0 and max_abs <= 1e-3:
        print("\nGATE PASSED: 100% argmax + max-abs <= 1e-3")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
