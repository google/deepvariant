// Phase 5.5c PASS-flip diagnostic: extract the pileup image at a
// specific (chrom, pos, ref, alt) from an examples TFRecord and write
// it as a single (1, H, W, C) NHWC FP32 .npy, where H/W/C are read from
// the example's image/shape feature (WGS 100x221x7 when absent). Pixel
// encoding matches call_variants ((src - 128) / 128).
//
// Used to byte-compare our pileup image against Docker's at the same
// site, isolating "inference drift" from "different pileup-image
// inputs".
//
// Usage:
//   extract_pileup_at_pos <examples.tfrecord[@N]> <out.npy> <chrom> <start_1based> <ref> <alt>
//
// Notes:
//   start_1based is the conventional VCF coordinate (1-based);
//   internally we compare against variant.start() which is 0-based.

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "deepvariant/native/tfrecord.h"

namespace {

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
    if ((tag & 7) != 2) break;
    uint64_t seg_len = ReadVarint(buf, len, i);
    if (i + seg_len > len) break;
    if (field == 1) {
      const uint8_t* inner = buf + i;
      size_t j = 0;
      while (j < seg_len) {
        uint64_t itag = ReadVarint(inner, seg_len, j);
        if ((itag & 7) != 2) break;
        uint64_t ilen = ReadVarint(inner, seg_len, j);
        if (j + ilen > seg_len) break;
        if ((itag >> 3) == 1) {
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

// Parse top-level tf.train.Example, return image bytes, variant bytes, shape.
struct ExampleParts {
  std::string image_encoded;
  std::string variant_encoded;
  std::vector<int64_t> image_shape;  // [H, W, C] when present & well-formed.
  bool image_shape_present = false;  // true if the feature key was seen at all,
                                     // independent of whether it decoded to 3.
};

ExampleParts ParseExample(const std::string& payload) {
  ExampleParts out;
  const uint8_t* buf = reinterpret_cast<const uint8_t*>(payload.data());
  size_t n = payload.size(), i = 0;
  while (i < n) {
    uint64_t tag = ReadVarint(buf, n, i);
    if ((tag & 7) != 2) break;
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
      std::string key, value_bytes;
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
      } else if (key == "variant/encoded") {
        out.variant_encoded = ExtractBytesListFirst(
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

// Variant proto field numbers (third_party/nucleus/protos/variants.proto):
//   reference_name = 14 (string), start = 16 (int64), end = 13 (int64),
//   reference_bases = 6 (string), alternate_bases = 7 (repeated string).
bool VariantMatches(const std::string& payload, const std::string& want_chrom,
                    int64_t want_start_0b, const std::string& want_ref,
                    const std::string& want_alt) {
  const uint8_t* buf = reinterpret_cast<const uint8_t*>(payload.data());
  size_t n = payload.size(), i = 0;
  std::string chrom, ref;
  std::vector<std::string> alts;
  int64_t start = -1;
  while (i < n) {
    uint64_t tag = ReadVarint(buf, n, i);
    uint32_t field = static_cast<uint32_t>(tag >> 3);
    uint32_t wire = static_cast<uint32_t>(tag & 7);
    if (wire == 0) {
      uint64_t v = ReadVarint(buf, n, i);
      if (field == 16) start = static_cast<int64_t>(v);
    } else if (wire == 2) {
      uint64_t seg_len = ReadVarint(buf, n, i);
      if (i + seg_len > n) break;
      if (field == 14) {
        chrom.assign(reinterpret_cast<const char*>(buf + i), seg_len);
      } else if (field == 6) {
        ref.assign(reinterpret_cast<const char*>(buf + i), seg_len);
      } else if (field == 7) {
        alts.emplace_back(reinterpret_cast<const char*>(buf + i), seg_len);
      }
      i += seg_len;
    } else if (wire == 5) {
      i += 4;
    } else if (wire == 1) {
      i += 8;
    } else {
      break;
    }
  }
  if (chrom != want_chrom) return false;
  if (start != want_start_0b) return false;
  if (ref != want_ref) return false;
  for (const auto& a : alts) if (a == want_alt) return true;
  return false;
}

bool WriteNpyFp32(const std::string& path, const float* data,
                   int N, int H, int W, int C) {
  std::ofstream f(path, std::ios::binary);
  if (!f) return false;
  std::string header =
      "{'descr': '<f4', 'fortran_order': False, 'shape': (" +
      std::to_string(N) + ", " + std::to_string(H) + ", " +
      std::to_string(W) + ", " + std::to_string(C) + "), }";
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
  size_t n_bytes = (size_t)N * H * W * C * sizeof(float);
  f.write(reinterpret_cast<const char*>(data), n_bytes);
  return f.good();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 7) {
    std::fprintf(stderr,
        "usage: %s <examples.tfrecord[@N]> <out.npy> <chrom> <start_1based> "
        "<ref> <alt>\n", argv[0]);
    return 2;
  }
  const std::string tfr_path = argv[1];
  const std::string out_path = argv[2];
  const std::string chrom = argv[3];
  // Checked parse of the 1-based start position. An unchecked strtoll() would
  // silently turn bad input ("abc", "", "12x") into position 0, scanning the
  // whole TFRecord and matching nothing — a confusing no-op. Reject anything
  // that isn't a complete, in-range positive integer.
  errno = 0;
  char* end = nullptr;
  const long long start_1based = std::strtoll(argv[4], &end, 10);
  if (errno != 0 || end == argv[4] || *end != '\0' || start_1based < 1) {
    std::fprintf(stderr, "error: <start_1based> must be a positive integer, "
                 "got \"%s\"\n", argv[4]);
    return 2;
  }
  const int64_t start_0b = static_cast<int64_t>(start_1based) - 1;
  const std::string ref = argv[5];
  const std::string alt = argv[6];

  auto reader = deepvariant::TFRecordReader::New(tfr_path);
  if (!reader) {
    std::fprintf(stderr, "cannot open %s\n", tfr_path.c_str());
    return 1;
  }

  int found = 0;
  long scanned = 0;
  while (reader->GetNext()) {
    ++scanned;
    auto p = ParseExample(reader->record());
    if (!VariantMatches(p.variant_encoded, chrom, start_0b, ref, alt)) continue;
    if (p.image_encoded.empty()) continue;
    // Geometry comes from the matched example's image/shape so non-WGS models
    // (WES/PacBio/ONT) work; fall back to WGS only when the feature is absent.
    // A present-but-non-3-D or out-of-range shape is an error, not a WGS guess.
    int H = 100, W = 221, C = 7;
    if (p.image_shape_present) {
      if (p.image_shape.size() != 3) {
        std::fprintf(stderr,
            "record %ld: image/shape has %zu values, expected 3 (H,W,C)\n",
            scanned, p.image_shape.size());
        return 1;
      }
      constexpr int64_t kMaxDim = 100000;
      const int64_t h = p.image_shape[0], w = p.image_shape[1],
                    c = p.image_shape[2];
      if (h <= 0 || w <= 0 || c <= 0 ||
          h > kMaxDim || w > kMaxDim || c > kMaxDim) {
        std::fprintf(stderr,
            "record %ld: image/shape %lldx%lldx%lld out of range (1..%lld)\n",
            scanned, (long long)h, (long long)w, (long long)c,
            (long long)kMaxDim);
        return 1;
      }
      H = static_cast<int>(h);
      W = static_cast<int>(w);
      C = static_cast<int>(c);
    }
    const int64_t kElem = static_cast<int64_t>(H) * W * C;
    std::vector<float> img(kElem);
    if ((int64_t)p.image_encoded.size() == kElem) {
      const uint8_t* src = reinterpret_cast<const uint8_t*>(p.image_encoded.data());
      constexpr float inv = 1.0f / 128.0f;
      for (int64_t j = 0; j < kElem; ++j) {
        img[j] = (static_cast<float>(src[j]) - 128.0f) * inv;
      }
    } else if ((int64_t)p.image_encoded.size() == kElem * 4) {
      std::memcpy(img.data(), p.image_encoded.data(),
                   (size_t)kElem * sizeof(float));
    } else {
      std::fprintf(stderr, "record %ld: bad image size %zu (expected %lld or "
                   "%lld for %dx%dx%d)\n", scanned, p.image_encoded.size(),
                   (long long)kElem, (long long)kElem * 4, H, W, C);
      continue;
    }
    if (!WriteNpyFp32(out_path, img.data(), 1, H, W, C)) {
      std::fprintf(stderr, "write failed: %s\n", out_path.c_str());
      return 1;
    }
    std::printf("MATCH at record %ld → wrote %s\n", scanned, out_path.c_str());
    ++found;
    break;
  }
  if (found == 0) {
    std::fprintf(stderr,
        "no match for %s:%lld %s>%s after %ld records\n",
        chrom.c_str(), (long long)start_0b + 1,
        ref.c_str(), alt.c_str(), scanned);
    return 1;
  }
  return 0;
}
