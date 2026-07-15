// Phase 5.5d — libstdc++-compatible std::shuffle for std::mt19937_64.
//
// Background: std::shuffle is implementation-defined; libc++ (Apple
// Clang) and libstdc++ (GCC, Docker) produce DIFFERENT sequences for
// the same input + generator state. This shows up as different
// pileup-image read selection in make_examples → different model input
// → 1.13 % FILTER drift vs `google/deepvariant:1.10.0` Docker on chr20.
//
// The cause is twofold:
//   1. Different Fisher–Yates iteration direction (forward vs
//      backward), hence different uniform_int call sequences.
//   2. Different `uniform_int_distribution<uint64_t>` algorithms —
//      libstdc++ 12 uses Lemire's nearly-divisionless method with
//      128-bit math; libc++ uses a rejection-sampling cousin.
//
// `LibstdcxxShuffle` reproduces libstdc++ 12's std::shuffle bit-for-bit
// for `std::vector<T>` with a `std::mt19937_64` generator (verified
// against `gcc:12` Docker on a 203-element vector with seed 2101079370 —
// first 20 + last 5 indices match exactly).
//
// Reference: libstdc++-v3/include/bits/stl_algo.h `shuffle`,
//            libstdc++-v3/include/bits/uniform_int_dist.h `_S_nd`.

#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace deepvariant {
namespace dv_shuffle {

// Lemire's nearly-divisionless uniform [0, range) using a 64-bit URBG.
// Mirrors libstdc++ uniform_int_distribution::_S_nd<__int128, …>(g, range).
inline uint64_t LemireUniformU64(std::mt19937_64& g, uint64_t range) {
  __extension__ typedef unsigned __int128 u128;
  u128 product = (u128)g() * (u128)range;
  uint64_t low = (uint64_t)product;
  if (low < range) {
    uint64_t threshold = (-range) % range;
    while (low < threshold) {
      product = (u128)g() * (u128)range;
      low = (uint64_t)product;
    }
  }
  return (uint64_t)(product >> 64);
}

// __gen_two_uniform_ints from stl_algo.h:
//   x = uniform_int(0, b0*b1 - 1)(g) → pair (x / b1, x % b1)
inline std::pair<uint64_t, uint64_t>
GenTwoUniformInts(std::mt19937_64& g, uint64_t b0, uint64_t b1) {
  uint64_t x = LemireUniformU64(g, b0 * b1);
  return {x / b1, x % b1};
}

// libstdc++ 12 std::shuffle — fast path for mt19937_64 (always taken
// when urange² < UINT64_MAX, i.e. always for our pileup sizes).
template <typename It>
inline void Shuffle(It first, It last, std::mt19937_64& g) {
  using DistanceType = typename std::iterator_traits<It>::difference_type;
  const DistanceType n = last - first;
  if (n < 2) return;
  uint64_t urange = (uint64_t)n;
  It i = first + 1;
  // If urange is even, swap count is uneven → handle leading swap solo.
  if ((urange % 2) == 0) {
    uint64_t r = LemireUniformU64(g, 2);     // uniform_int(0, 1)
    std::iter_swap(i, first + (DistanceType)r);
    ++i;
  }
  while (i != last) {
    const uint64_t swap_range = (uint64_t)(i - first) + 1;
    auto pp = GenTwoUniformInts(g, swap_range, swap_range + 1);
    std::iter_swap(i, first + (DistanceType)pp.first);
    ++i;
    std::iter_swap(i, first + (DistanceType)pp.second);
    ++i;
  }
}

}  // namespace dv_shuffle
}  // namespace deepvariant
