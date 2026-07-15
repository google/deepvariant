// Phase 5.5d/3 — NumPy-compatible MT19937 + bounded_lemire_uint32 +
// Algorithm-R reservoir sampling, ported from numpy/random/src/mt19937/
// (NumPy 1.24, the version inside `google/deepvariant:1.10.0` Docker).
//
// Used by `make_examples_main.cc` to subsample reads per partition
// (max_reads_per_partition=1500) using the same RNG sequence as
// `np.random.RandomState(seed).randint(0, n)` and the same Algorithm R
// reservoir as `numpy/utils.py::reservoir_sample`. This is what closes
// the chr20 DP mismatch at high-coverage outlier sites — same input set
// reaches AlleleCounter on both sides.
//
// Verified against `np.random.RandomState(2101079370)` test vectors:
//   randint(0, 1000) ×10  = 940, 785, 301, 77, 558, 250, 667, 359, 899, 910
//   randint(0, i+1) i=0..19 = 0, 0, 1, 1, 2, 3, 3, 6, 7, 5,
//                              10, 7, 5, 5, 9, 7, 2, 3, 9, 9

#pragma once

#include <cstdint>
#include <vector>

namespace deepvariant {
namespace npr {

// Standard MT19937 (Matsumoto-Nishimura 1998) — same engine NumPy uses
// for legacy `RandomState`. State = 624 × uint32. Tempering output.
class NumpyMt19937 {
 public:
  static constexpr int kStateLen = 624;
  static constexpr int kMid = 397;
  static constexpr uint32_t kMatrixA = 0x9908b0dfUL;
  static constexpr uint32_t kUpperMask = 0x80000000UL;
  static constexpr uint32_t kLowerMask = 0x7fffffffUL;

  // Seed via the canonical `init_genrand` (a.k.a. mt19937_seed). The
  // 1812433253 multiplier is Matsumoto-Nishimura's; NumPy uses the same.
  explicit NumpyMt19937(uint32_t seed) {
    state_[0] = seed;
    for (int i = 1; i < kStateLen; ++i) {
      state_[i] =
          (1812433253UL * (state_[i - 1] ^ (state_[i - 1] >> 30)) + i);
    }
    pos_ = kStateLen;
  }

  uint32_t NextUint32() {
    if (pos_ >= kStateLen) Generate();
    uint32_t y = state_[pos_++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
  }

 private:
  void Generate() {
    static constexpr uint32_t mag01[2] = {0, kMatrixA};
    int i;
    for (i = 0; i < kStateLen - kMid; ++i) {
      uint32_t y = (state_[i] & kUpperMask) | (state_[i + 1] & kLowerMask);
      state_[i] = state_[i + kMid] ^ (y >> 1) ^ mag01[y & 1];
    }
    for (; i < kStateLen - 1; ++i) {
      uint32_t y = (state_[i] & kUpperMask) | (state_[i + 1] & kLowerMask);
      state_[i] = state_[i + (kMid - kStateLen)] ^ (y >> 1) ^ mag01[y & 1];
    }
    uint32_t y =
        (state_[kStateLen - 1] & kUpperMask) | (state_[0] & kLowerMask);
    state_[kStateLen - 1] = state_[kMid - 1] ^ (y >> 1) ^ mag01[y & 1];
    pos_ = 0;
  }

  uint32_t state_[kStateLen];
  int pos_;
};

// NumPy `random_interval(bg, max)` — uniform integer in [0, max] inclusive.
// Mirrors `numpy/random/src/distributions/distributions.c::random_interval`:
// build the next-power-of-2 mask ≥ max, draw a u32, mask it, accept iff
// ≤ max. NOT Lemire — that's used elsewhere in NumPy (e.g.,
// `Generator.integers`), but the legacy `RandomState.randint` path goes
// through `random_interval`.
inline uint32_t NumpyRandomIntervalU32(NumpyMt19937& g, uint32_t max_inc) {
  if (max_inc == 0) return 0;
  uint32_t mask = max_inc;
  mask |= mask >> 1;
  mask |= mask >> 2;
  mask |= mask >> 4;
  mask |= mask >> 8;
  mask |= mask >> 16;
  uint32_t value;
  do {
    value = g.NextUint32() & mask;
  } while (value > max_inc);
  return value;
}

// `np.random.RandomState(seed).randint(0, n)` — returns uniform [0, n).
inline uint32_t RandintU32(NumpyMt19937& g, uint32_t n) {
  if (n == 0) return 0;
  return NumpyRandomIntervalU32(g, n - 1);
}

// Algorithm R reservoir sampling, mirror of
// `third_party/nucleus/util/utils.py::reservoir_sample`:
//
//   sample = []
//   for i, item in enumerate(iterable):
//       if len(sample) < k:
//           sample.append(item)
//       else:
//           j = random.randint(0, i + 1)   # uniform [0, i]
//           if j < k:
//               sample[j] = item
//   return sample
//
// `k` is the cap; `iterable` is anything with stable iteration order
// (we keep a vector of pointers to avoid copying T). Returns the
// retained pointers in the order they sit in the reservoir at the end
// — same as upstream.
template <typename T>
std::vector<const T*> ReservoirSamplePtrs(
    const std::vector<T>& items, size_t k, NumpyMt19937& gen) {
  std::vector<const T*> sample;
  sample.reserve(std::min(items.size(), k));
  for (size_t i = 0; i < items.size(); ++i) {
    if (sample.size() < k) {
      sample.push_back(&items[i]);
    } else {
      // randint(0, i + 1) — uniform [0, i] inclusive.
      uint32_t j = RandintU32(gen, (uint32_t)(i + 1));
      if (j < k) sample[j] = &items[i];
    }
  }
  return sample;
}

}  // namespace npr
}  // namespace deepvariant
