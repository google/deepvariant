// Profiling tool: extract the first N pileup images from a TFRecord (or
// `name@N` shard spec) and write them as a NumPy `.npy` array of shape
// (N, H, W, C) FP32 NHWC, where H/W/C are read from the example's
// image/shape feature (WGS 100x221x7 when absent).  Pixel encoding mirrors
// call_variants:
//   uint8 src → (src - 128) / 128.0 → FP32
// or a passthrough when the input is already FP32.
//
// Used by Phase 5.5c per-layer drift profiling: produces a real-data
// `_input.npy` that `dump_tf_per_layer.py` (Docker) and
// `debug_metal --compare-to-reference` both consume.
//
// usage:
//   extract_pileup_npy <examples.tfrecord[@N]> <out.npy> [count=64]

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "deepvariant/native/tfrecord.h"

namespace {

// Minimal protobuf wire decoders — ported from call_variants_main.cc's
// anonymous namespace (we don't link the runtime here).

uint64_t ReadVarint(const uint8_t* buf, size_t len, size_t& i) {
  uint64_t val = 0;
  int shift = 0;
  while (i < len) {
    uint8_t b = buf[i++];
    val |= static_cast<uint64_t>(b & 0x7F) << shift;
    if (!(b & 0x80)) return val;
    shift += 7;
  }
  return val;
}

std::string ExtractBytesListFirst(const uint8_t* buf, size_t len) {
  size_t i = 0;
  while (i < len) {
    uint64_t tag = ReadVarint(buf, len, i);
    uint32_t field = static_cast<uint32_t>(tag >> 3);
    uint32_t wire  = static_cast<uint32_t>(tag & 7);
    if (wire != 2) break;
    uint64_t seg_len = ReadVarint(buf, len, i);
    if (i + seg_len > len) break;
    if (field == 1) {
      const uint8_t* inner = buf + i;
      size_t j = 0;
      while (j < seg_len) {
        uint64_t itag = ReadVarint(inner, seg_len, j);
        uint32_t ifield = static_cast<uint32_t>(itag >> 3);
        uint32_t iwire  = static_cast<uint32_t>(itag & 7);
        if (iwire != 2) break;
        uint64_t ilen = ReadVarint(inner, seg_len, j);
        if (j + ilen > seg_len) break;
        if (ifield == 1) {
          return std::string(reinterpret_cast<const char*>(inner + j), ilen);
        }
        j += ilen;
      }
      return {};
    }
    i += seg_len;
  }
  return {};
}

// Decode a tf.train.Feature message holding an Int64List into its values.
// Handles both packed (proto3 default) and unpacked repeated-int64 encodings.
// Returns empty if the Feature is not an Int64List.
std::vector<int64_t> ParseInt64List(const uint8_t* buf, size_t len) {
  std::vector<int64_t> out;
  size_t i = 0;
  while (i < len) {
    uint64_t tag = ReadVarint(buf, len, i);
    uint32_t field = static_cast<uint32_t>(tag >> 3);
    uint32_t wire  = static_cast<uint32_t>(tag & 7);
    if (wire == 2) {
      uint64_t seg_len = ReadVarint(buf, len, i);
      if (i + seg_len > len) break;
      if (field == 3) {  // Feature.int64_list
        const uint8_t* sub = buf + i;
        size_t si = 0;
        while (si < seg_len) {
          uint64_t vtag = ReadVarint(sub, seg_len, si);
          uint32_t vfield = static_cast<uint32_t>(vtag >> 3);
          uint32_t vwire  = static_cast<uint32_t>(vtag & 7);
          if (vfield == 1 && vwire == 2) {        // packed values
            uint64_t plen = ReadVarint(sub, seg_len, si);
            if (si + plen > seg_len) break;
            size_t pend = si + plen;
            while (si < pend) {
              // Bound reads by pend (the packed-blob end), not seg_len, so a
              // truncated trailing varint can't run into the rest of the message.
              out.push_back(static_cast<int64_t>(ReadVarint(sub, pend, si)));
            }
          } else if (vfield == 1 && vwire == 0) {  // single unpacked value
            out.push_back(static_cast<int64_t>(ReadVarint(sub, seg_len, si)));
          } else {
            break;
          }
        }
        return out;
      }
      i += seg_len;
    } else if (wire == 0) {
      ReadVarint(buf, len, i);
    } else if (wire == 5) {
      i += 4;
    } else if (wire == 1) {
      i += 8;
    } else {
      break;
    }
  }
  return out;
}

struct ParsedExample {
  std::string image_encoded;
  std::vector<int64_t> image_shape;  // [H, W, C] when present & well-formed.
  bool image_shape_present = false;  // true if the feature key was seen at all,
                                     // independent of whether it decoded to 3.
};

ParsedExample ParseExample(const std::string& payload) {
  ParsedExample out;
  const uint8_t* buf = reinterpret_cast<const uint8_t*>(payload.data());
  size_t n = payload.size();
  size_t i = 0;
  while (i < n) {
    uint64_t tag = ReadVarint(buf, n, i);
    uint32_t wire = tag & 7;
    if (wire != 2) break;
    uint64_t seg_len = ReadVarint(buf, n, i);
    if (i + seg_len > n) break;
    const uint8_t* feat_buf = buf + i;
    size_t feat_len = seg_len;
    i += seg_len;
    size_t fi = 0;
    while (fi < feat_len) {
      uint64_t ftag = ReadVarint(feat_buf, feat_len, fi);
      if ((ftag & 7) != 2) break;
      uint64_t entry_len = ReadVarint(feat_buf, feat_len, fi);
      if (fi + entry_len > feat_len) break;
      const uint8_t* entry = feat_buf + fi;
      fi += entry_len;
      std::string key;
      std::string value_bytes;
      size_t ei = 0;
      while (ei < entry_len) {
        uint64_t etag = ReadVarint(entry, entry_len, ei);
        uint32_t efd = etag >> 3;
        if ((etag & 7) != 2) break;
        uint64_t elen = ReadVarint(entry, entry_len, ei);
        if (ei + elen > entry_len) break;
        if (efd == 1) {
          key.assign(reinterpret_cast<const char*>(entry + ei), elen);
        } else if (efd == 2) {
          value_bytes.assign(reinterpret_cast<const char*>(entry + ei), elen);
        }
        ei += elen;
      }
      if (key == "image/encoded" || key == "image") {
        out.image_encoded = ExtractBytesListFirst(
            reinterpret_cast<const uint8_t*>(value_bytes.data()),
            value_bytes.size());
      } else if (key == "image/shape") {
        out.image_shape_present = true;
        out.image_shape = ParseInt64List(
            reinterpret_cast<const uint8_t*>(value_bytes.data()),
            value_bytes.size());
      }
    }
  }
  return out;
}

// Write a (N, H, W, C) FP32 NHWC array to NumPy v1 .npy.
bool WriteNpyFp32NHWC(const std::string& path, int N, int H, int W, int C,
                      const float* data) {
  std::ofstream f(path, std::ios::binary);
  if (!f) return false;

  std::string header =
      "{'descr': '<f4', 'fortran_order': False, 'shape': (" +
      std::to_string(N) + ", " + std::to_string(H) + ", " +
      std::to_string(W) + ", " + std::to_string(C) + "), }";
  // Pad header so total prefix (10 bytes magic+version+len + header + 1 \n)
  // is a multiple of 64 — required by the .npy format.
  while (((10 + header.size() + 1) % 64) != 0) header.push_back(' ');
  header.push_back('\n');

  const char magic[6] = {'\x93','N','U','M','P','Y'};
  f.write(magic, 6);
  uint8_t major = 1, minor = 0;
  f.write(reinterpret_cast<const char*>(&major), 1);
  f.write(reinterpret_cast<const char*>(&minor), 1);
  uint16_t hl = static_cast<uint16_t>(header.size());
  f.write(reinterpret_cast<const char*>(&hl), 2);
  f.write(header.data(), header.size());

  const size_t n_bytes =
      static_cast<size_t>(N) * H * W * C * sizeof(float);
  f.write(reinterpret_cast<const char*>(data), n_bytes);
  return f.good();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 4) {
    std::fprintf(stderr,
        "usage: %s <examples.tfrecord[@N]> <out.npy> [count=64]\n", argv[0]);
    return 2;
  }
  const std::string tfr_path = argv[1];
  const std::string out_path = argv[2];
  const int count = (argc >= 4) ? std::atoi(argv[3]) : 64;
  if (count <= 0 || count > 100000) {
    std::fprintf(stderr, "bad count=%d\n", count);
    return 2;
  }

  auto reader = deepvariant::TFRecordReader::New(tfr_path);
  if (!reader) {
    std::fprintf(stderr, "cannot open %s\n", tfr_path.c_str());
    return 1;
  }

  // Geometry is taken from the first record's image/shape feature so the tool
  // adapts to any model type (WES/PacBio/ONT differ from WGS). The whole batch
  // is packed into one (N, H, W, C) array, so later records must share the
  // geometry — a mismatch is an error rather than a silent overwrite.
  int H = 0, W = 0, C = 0;
  int64_t kElemPerImg = 0;
  std::vector<float> all;
  int n_loaded = 0;
  for (int i = 0; i < count; ++i) {
    if (!reader->GetNext()) {
      std::fprintf(stderr, "EOF after %d records\n", i);
      break;
    }
    const ParsedExample ex = ParseExample(reader->record());
    if (i == 0) {
      if (!ex.image_shape_present) {
        H = 100; W = 221; C = 7;  // WGS fallback only when the feature is absent.
        std::fprintf(stderr,
            "warning: record 0 has no image/shape; assuming WGS %dx%dx%d\n",
            H, W, C);
      } else if (ex.image_shape.size() != 3) {
        // Present but not a 3-D [H,W,C] shape — don't silently guess WGS.
        std::fprintf(stderr,
            "record 0: image/shape has %zu values, expected 3 (H,W,C)\n",
            ex.image_shape.size());
        return 1;
      } else {
        // Validate each int64 dim is a sane positive value before narrowing to
        // int, so a corrupt shape can't wrap to garbage or request an absurd
        // allocation.
        constexpr int64_t kMaxDim = 100000;
        const int64_t h = ex.image_shape[0], w = ex.image_shape[1],
                      c = ex.image_shape[2];
        if (h <= 0 || w <= 0 || c <= 0 ||
            h > kMaxDim || w > kMaxDim || c > kMaxDim) {
          std::fprintf(stderr,
              "record 0: image/shape %lldx%lldx%lld out of range (1..%lld)\n",
              (long long)h, (long long)w, (long long)c, (long long)kMaxDim);
          return 1;
        }
        H = static_cast<int>(h);
        W = static_cast<int>(w);
        C = static_cast<int>(c);
      }
      kElemPerImg = static_cast<int64_t>(H) * W * C;
      // count is already bounded to (0, 100000]; guard the batch allocation.
      if (kElemPerImg <= 0 || kElemPerImg > INT64_MAX / count) {
        std::fprintf(stderr, "record 0: batch too large (%lld elems x %d)\n",
                     (long long)kElemPerImg, count);
        return 1;
      }
      all.resize(static_cast<size_t>(count) * kElemPerImg);
    } else if (ex.image_shape.size() == 3 &&
               (ex.image_shape[0] != H || ex.image_shape[1] != W ||
                ex.image_shape[2] != C)) {
      std::fprintf(stderr,
          "record %d: image/shape %lldx%lldx%lld differs from batch %dx%dx%d; "
          "cannot pack heterogeneous geometry into one array\n",
          i, static_cast<long long>(ex.image_shape[0]),
          static_cast<long long>(ex.image_shape[1]),
          static_cast<long long>(ex.image_shape[2]), H, W, C);
      return 1;
    }
    const std::string& img = ex.image_encoded;
    float* dst = all.data() + static_cast<size_t>(i) * kElemPerImg;
    if (static_cast<int64_t>(img.size()) == kElemPerImg) {
      // uint8 → (x - 128) / 128 — same path as call_variants.
      const uint8_t* src = reinterpret_cast<const uint8_t*>(img.data());
      constexpr float inv = 1.0f / 128.0f;
      for (int64_t j = 0; j < kElemPerImg; ++j) {
        dst[j] = (static_cast<float>(src[j]) - 128.0f) * inv;
      }
    } else if (static_cast<int64_t>(img.size()) == kElemPerImg * 4) {
      std::memcpy(dst, img.data(),
                  static_cast<size_t>(kElemPerImg) * sizeof(float));
    } else {
      std::fprintf(stderr,
          "record %d: bad image size %zu (expected %lld or %lld)\n",
          i, img.size(),
          static_cast<long long>(kElemPerImg),
          static_cast<long long>(kElemPerImg * 4));
      return 1;
    }
    ++n_loaded;
  }

  if (n_loaded == 0) {
    // Empty input / immediate EOF: don't emit a bogus (0,0,0,0) .npy.
    std::fprintf(stderr, "no records read from %s; nothing written\n",
                 tfr_path.c_str());
    return 1;
  }

  if (!WriteNpyFp32NHWC(out_path, n_loaded, H, W, C, all.data())) {
    std::fprintf(stderr, "failed to write %s\n", out_path.c_str());
    return 1;
  }
  std::printf("wrote %d images to %s (shape %d×%d×%d×%d, %.1f MB)\n",
              n_loaded, out_path.c_str(), n_loaded, H, W, C,
              n_loaded * kElemPerImg * 4.0 / (1024.0 * 1024.0));
  return 0;
}
