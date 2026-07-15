// Loader for the `.dvw` weight bundle format produced by
// `tools/conversion/extract_weights.py`. mmap-backed, zero-copy access
// to FP32 tensors keyed by name.
//
// Used by the Phase 5.5 Metal/BNNS inference path to load model weights
// at runtime without depending on TensorFlow, coremltools, or any proto
// runtime.
//
// File layout (all integers little-endian, see extract_weights.py):
//
//     magic[4] = 'DVW1'
//     version[4] = 1
//     n_tensors[4]
//     for each tensor (sorted by name for determinism):
//         name_len[4]
//         name[name_len]    (utf-8)
//         dtype[1]          (1 = DT_FLOAT)
//         ndim[1]
//         shape[ndim*4]     (uint32 le)
//         offset[8]         (into payload)
//         n_bytes[8]
//     payload: concatenated raw FP32 LE bytes
//
// Threadsafe for read-only access after Open().
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace deepvariant {

struct DvwTensor {
  // Raw pointer into the mmap'd file. Valid as long as the parent
  // DvwWeights object is alive.
  const float* data = nullptr;
  std::vector<uint32_t> shape;
  size_t n_elements = 0;  // product(shape)
  size_t n_bytes = 0;     // n_elements * sizeof(float)
};

class DvwWeights {
 public:
  // Open and parse a .dvw file. Returns nullptr on any error
  // (file missing, bad magic, truncated table, etc.).
  static std::unique_ptr<DvwWeights> Open(const std::string& path);

  ~DvwWeights();

  // Look up a tensor by its source name (e.g.
  // "layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE").
  // Returns nullptr if absent. The returned pointer is owned by `this`.
  const DvwTensor* Get(const std::string& name) const;

  // Iterate all tensor names (sorted as on-disk order).
  const std::vector<std::string>& Names() const { return names_; }

  uint32_t Version() const { return version_; }

  DvwWeights(const DvwWeights&) = delete;
  DvwWeights& operator=(const DvwWeights&) = delete;

 private:
  DvwWeights();

  // Owned mmap mapping.
  void* map_addr_ = nullptr;
  size_t map_size_ = 0;

  uint32_t version_ = 0;
  std::vector<std::string> names_;
  std::unordered_map<std::string, DvwTensor> by_name_;
};

}  // namespace deepvariant
