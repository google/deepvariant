"""TF-free TensorBundle reader.

Reads `variables/variables.index` (a leveldb-style sstable) and the
matching `variables/variables.data-NNNNN-of-MMMMM` shards to recover
the tensor weights of a TensorFlow SavedModel — without importing TF.

Format reference: tensorflow/core/util/tensor_bundle/tensor_bundle.h
                  tensorflow/core/lib/io/{table.cc, format.cc, block.cc}
                  github.com/google/leveldb/doc/table_format.md

High-level flow:
    bundle = read_bundle(prefix="variables/variables")
    for name, entry in bundle.entries.items():
        # entry.dtype, entry.shape, entry.shard_id, entry.offset, entry.size
        arr = bundle.read_tensor(name)   # numpy.ndarray
"""

from __future__ import annotations

import struct
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Bring the generated protobuf bindings onto sys.path.
_GEN = Path(__file__).resolve().parent / "Generated"
if str(_GEN) not in sys.path:
    sys.path.insert(0, str(_GEN))


# ---------------------------------------------------------------------------
# Constants from leveldb / TF table format
# ---------------------------------------------------------------------------

# tensorflow/core/lib/io/format.cc: kTableMagicNumber
_TABLE_MAGIC = 0xDB4775248B80FB57
_FOOTER_SIZE = 48  # legacy TF/leveldb footer: 2*BlockHandle (≤20B each) + 8B magic
_BLOCK_TRAILER_SIZE = 5  # 1 byte type + 4 bytes CRC

_NO_COMPRESSION = 0
_SNAPPY_COMPRESSION = 1


# ---------------------------------------------------------------------------
# Varints
# ---------------------------------------------------------------------------

def _read_varint(buf: bytes, i: int) -> tuple[int, int]:
    """Read a base-128 varint starting at offset i. Returns (value, new_i)."""
    val = 0
    shift = 0
    while True:
        b = buf[i]
        i += 1
        val |= (b & 0x7F) << shift
        if not (b & 0x80):
            return val, i
        shift += 7
        if shift > 63:
            raise RuntimeError("varint too long")


# ---------------------------------------------------------------------------
# Footer + BlockHandle
# ---------------------------------------------------------------------------

@dataclass
class BlockHandle:
    offset: int
    size: int


def _decode_block_handle(buf: bytes, i: int) -> tuple[BlockHandle, int]:
    offset, i = _read_varint(buf, i)
    size, i = _read_varint(buf, i)
    return BlockHandle(offset, size), i


def _read_footer(path: Path) -> tuple[BlockHandle, BlockHandle]:
    """Read the last 48 bytes of `path` and decode metaindex + index handles."""
    with open(path, "rb") as f:
        f.seek(-_FOOTER_SIZE, 2)
        footer = f.read(_FOOTER_SIZE)
    if len(footer) != _FOOTER_SIZE:
        raise RuntimeError(f"truncated footer in {path}")
    magic = struct.unpack("<Q", footer[-8:])[0]
    if magic != _TABLE_MAGIC:
        raise RuntimeError(
            f"bad table magic in {path}: "
            f"got 0x{magic:016x}, want 0x{_TABLE_MAGIC:016x}"
        )
    metaindex, i = _decode_block_handle(footer, 0)
    index, _ = _decode_block_handle(footer, i)
    return metaindex, index


# ---------------------------------------------------------------------------
# Block reader
# ---------------------------------------------------------------------------

def _snappy_decompress(data: bytes) -> bytes:
    """Pure-Python Snappy decompressor.

    Format: varint uncompressed_length, then a stream of tag-prefixed
    sub-blocks. Tag bits 0-1 = element type:
        0b00 LITERAL   ((tag>>2)+1 bytes; or extended for ≥60)
        0b01 COPY_1    1-byte offset (high 3 bits in tag)
        0b10 COPY_2    2-byte LE offset
        0b11 COPY_4    4-byte LE offset
    See google/snappy/blob/main/format_description.txt.
    """
    expected_len, i = _read_varint(data, 0)
    out = bytearray()
    n = len(data)
    while i < n:
        tag = data[i]
        i += 1
        kind = tag & 0x3
        if kind == 0:  # LITERAL
            ln = tag >> 2
            if ln < 60:
                ln += 1
            else:
                extra = ln - 59  # 1, 2, 3 or 4
                ln = (
                    int.from_bytes(data[i:i + extra], "little") + 1
                )
                i += extra
            out += data[i:i + ln]
            i += ln
        else:
            if kind == 1:  # COPY_1
                ln = ((tag >> 2) & 0x7) + 4
                offset = ((tag >> 5) << 8) | data[i]
                i += 1
            elif kind == 2:  # COPY_2
                ln = (tag >> 2) + 1
                offset = int.from_bytes(data[i:i + 2], "little")
                i += 2
            else:  # COPY_4
                ln = (tag >> 2) + 1
                offset = int.from_bytes(data[i:i + 4], "little")
                i += 4
            # Snappy may copy past current end (overlapping copy);
            # naive byte-by-byte handles that correctly.
            base = len(out) - offset
            for k in range(ln):
                out.append(out[base + k])
    if len(out) != expected_len:
        raise RuntimeError(
            f"snappy: decoded {len(out)} bytes, expected {expected_len}"
        )
    return bytes(out)


def _read_block(path: Path, handle: BlockHandle) -> bytes:
    """Read the data portion of a block (drops the 5-byte type+CRC trailer).

    No CRC verification (dev tool). Supports uncompressed + snappy —
    TF tensor bundles default to snappy.
    """
    with open(path, "rb") as f:
        f.seek(handle.offset)
        raw = f.read(handle.size + _BLOCK_TRAILER_SIZE)
    if len(raw) != handle.size + _BLOCK_TRAILER_SIZE:
        raise RuntimeError(
            f"truncated block at offset {handle.offset} in {path}"
        )
    body = raw[:handle.size]
    comp_type = raw[handle.size]
    if comp_type == _NO_COMPRESSION:
        return body
    if comp_type == _SNAPPY_COMPRESSION:
        return _snappy_decompress(body)
    raise RuntimeError(f"unknown block compression type {comp_type}")


def _iterate_block(block: bytes):
    """Yield (key, value) pairs from a leveldb-format data block.

    Block layout:
        records...
        restart_array: uint32 * num_restarts
        num_restarts:  uint32   (last 4 bytes)

    Each record:
        shared_key_len:    varint
        non_shared_key_len varint
        value_len:         varint
        non_shared_key:    bytes
        value:             bytes
    """
    if len(block) < 4:
        return
    num_restarts = struct.unpack("<I", block[-4:])[0]
    records_end = len(block) - 4 - 4 * num_restarts
    i = 0
    prev_key = b""
    while i < records_end:
        shared, i = _read_varint(block, i)
        non_shared, i = _read_varint(block, i)
        value_len, i = _read_varint(block, i)
        non_shared_key = block[i:i + non_shared]
        i += non_shared
        value = block[i:i + value_len]
        i += value_len
        key = prev_key[:shared] + non_shared_key
        yield key, value
        prev_key = key


# ---------------------------------------------------------------------------
# Bundle
# ---------------------------------------------------------------------------

# DataType -> (numpy dtype, byte size).
# Only the dtypes we actually expect in DV's Inception-v3 bundle.
_NUMPY_DTYPE = {
    1: (np.float32, 4),     # DT_FLOAT
    2: (np.float64, 8),     # DT_DOUBLE
    3: (np.int32, 4),       # DT_INT32
    9: (np.int64, 8),       # DT_INT64
    10: (np.bool_, 1),      # DT_BOOL
    14: (np.float32, 2),    # DT_BFLOAT16 (2 bytes on disk; widened to fp32 on read since bfloat16 isn't in numpy)
    17: (np.float16, 2),    # DT_HALF
}


@dataclass
class TensorEntry:
    name: str
    dtype: int
    shape: tuple[int, ...]
    shard_id: int
    offset: int
    size: int


class TensorBundle:
    """A read-only view of a TF tensor bundle (one .index + N .data shards).

    The .index file is a leveldb-format table:
        footer -> index_block -> data_block(s) -> records
    The records' keys are variable names (UTF-8) or empty (== header).
    Records' values are serialized BundleEntryProto / BundleHeaderProto.
    """

    def __init__(self, prefix: str | Path) -> None:
        self.prefix = Path(prefix)
        index_path = self.prefix.with_suffix(".index") if self.prefix.suffix \
            else Path(str(self.prefix) + ".index")
        if not index_path.exists():
            raise FileNotFoundError(f"{index_path} not found")
        self.index_path = index_path

        from tensorflow.core.protobuf import (  # type: ignore
            tensor_bundle_pb2,
        )
        self._tb_pb = tensor_bundle_pb2

        _, index_handle = _read_footer(index_path)
        index_block = _read_block(index_path, index_handle)

        self.header = None
        self.entries: dict[str, TensorEntry] = {}

        # The index block's records each map a "largest key in the data
        # block" -> BlockHandle (encoded as two varints) of that data block.
        # Iterate each data block and collect the actual records.
        for _largest_key, handle_bytes in _iterate_block(index_block):
            data_handle, _ = _decode_block_handle(handle_bytes, 0)
            data_block = _read_block(index_path, data_handle)
            for key, value in _iterate_block(data_block):
                if not key:
                    hdr = tensor_bundle_pb2.BundleHeaderProto()
                    hdr.ParseFromString(value)
                    self.header = hdr
                    continue
                entry = tensor_bundle_pb2.BundleEntryProto()
                entry.ParseFromString(value)
                name = key.decode("utf-8")
                shape = tuple(d.size for d in entry.shape.dim)
                self.entries[name] = TensorEntry(
                    name=name,
                    dtype=entry.dtype,
                    shape=shape,
                    shard_id=entry.shard_id,
                    offset=entry.offset,
                    size=entry.size,
                )

        if self.header is None:
            raise RuntimeError(
                f"{index_path}: no BundleHeaderProto entry (empty key)"
            )

    def shard_path(self, shard_id: int) -> Path:
        n = self.header.num_shards
        name = (
            f"{self.prefix.name}.data-{shard_id:05d}-of-{n:05d}"
        )
        return self.prefix.parent / name

    def read_raw(self, name: str) -> bytes:
        """Return the raw bytes of a bundle entry (no dtype interpretation)."""
        if name not in self.entries:
            raise KeyError(f"variable {name!r} not in bundle")
        e = self.entries[name]
        path = self.shard_path(e.shard_id)
        with open(path, "rb") as f:
            f.seek(e.offset)
            raw = f.read(e.size)
        if len(raw) != e.size:
            raise RuntimeError(
                f"truncated read of {name} from {path}: "
                f"want {e.size} bytes, got {len(raw)}"
            )
        return raw

    def read_tensor(self, name: str) -> np.ndarray:
        if name not in self.entries:
            raise KeyError(f"variable {name!r} not in bundle")
        e = self.entries[name]
        if e.dtype not in _NUMPY_DTYPE:
            raise NotImplementedError(
                f"dtype {e.dtype} not supported (variable {name})"
            )
        np_dtype, dtype_size = _NUMPY_DTYPE[e.dtype]
        path = self.shard_path(e.shard_id)
        with open(path, "rb") as f:
            f.seek(e.offset)
            raw = f.read(e.size)
        if len(raw) != e.size:
            raise RuntimeError(
                f"truncated read of {name} from {path}: "
                f"want {e.size} bytes, got {len(raw)}"
            )
        # Guard: the raw byte count must be a whole number of elements that
        # matches the declared shape. Otherwise np.frombuffer would silently
        # mis-decode (e.g. reinterpreting 2-byte values as 4-byte ones).
        expected_bytes = int(np.prod(e.shape)) * dtype_size
        if len(raw) != expected_bytes:
            raise ValueError(
                f"size mismatch for variable {name}: shape {e.shape} with "
                f"{dtype_size}-byte elements implies {expected_bytes} bytes, "
                f"but got {len(raw)} bytes"
            )
        if e.dtype == 14:  # DT_BFLOAT16: 2 raw bytes are the high 16 bits of fp32.
            u16 = np.frombuffer(raw, dtype=np.uint16)
            arr = (u16.astype(np.uint32) << 16).view(np.float32)
        else:
            arr = np.frombuffer(raw, dtype=np_dtype)
        if e.shape:
            arr = arr.reshape(e.shape)
        return arr

    def names(self) -> list[str]:
        return sorted(self.entries.keys())


def open_bundle(prefix: str | Path) -> TensorBundle:
    """Convenience entry — same as `TensorBundle(prefix)`."""
    return TensorBundle(prefix)


