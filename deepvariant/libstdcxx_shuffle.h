// Copyright 2026 Google LLC.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef LEARNING_GENOMICS_DEEPVARIANT_LIBSTDCXX_SHUFFLE_H_
#define LEARNING_GENOMICS_DEEPVARIANT_LIBSTDCXX_SHUFFLE_H_

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <random>
#include <utility>

namespace learning {
namespace genomics {
namespace deepvariant {
namespace dv_shuffle {

// Lemire's nearly-divisionless uniform [0, range) using a 64-bit URBG.
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

inline std::pair<uint64_t, uint64_t> GenTwoUniformInts(std::mt19937_64& g,
                                                       uint64_t b0,
                                                       uint64_t b1) {
  uint64_t x = LemireUniformU64(g, b0 * b1);
  return {x / b1, x % b1};
}

// libstdc++ 12 std::shuffle for mt19937_64.
template <typename It>
// NOLINTNEXTLINE(runtime/expensive_random_temporaries)
inline void Shuffle(It first, It last, std::mt19937_64& g) {
  using DistanceType = typename std::iterator_traits<It>::difference_type;
  const DistanceType n = last - first;
  if (n < 2) return;
  uint64_t urange = (uint64_t)n;
  It i = first + 1;
  if ((urange % 2) == 0) {
    uint64_t r = LemireUniformU64(g, 2);
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
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_LIBSTDCXX_SHUFFLE_H_
