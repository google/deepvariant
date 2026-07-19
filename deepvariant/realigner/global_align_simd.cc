/*
 * Copyright 2026 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

// global_align_simd.cc — SIMD-vectorized GlobalAlign driver
//
// This file provides the HWY_EXPORT dispatch for the SIMD DP and the
// striped-layout BackTrackBestAlignment implementation.

// clang-format off
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE \
  "deepvariant/realigner/global_align_simd.cc"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "hwy/aligned_allocator.h"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

#include "deepvariant/realigner/global_align_simd-inl.h"
#include "deepvariant/realigner/global_align_simd.h"
// clang-format on

// ============================================================================
// Non-SIMD code — compiled once.
// ============================================================================
#if HWY_ONCE

namespace learning {
namespace genomics {
namespace deepvariant {

// Export the SIMD functions for dynamic dispatch.
HWY_EXPORT(PopulateDpMatrixSimd);
HWY_EXPORT(SimdWidth);

namespace {

// Get the runtime SIMD width for int16.
size_t GetSimdWidth() {
  return HWY_DYNAMIC_DISPATCH(SimdWidth)();
}

// Access a value in the striped layout.
// i: 1-indexed query position (1..n)
// j: column index (0..m)
int16_t GetStriped(const int16_t* matrix, int i, int j, int S, size_t W,
                   size_t stripe_size) {
  int p = i - 1;  // 0-indexed query position
  int seg = p % S;
  size_t lane = static_cast<size_t>(p / S);
  return matrix[static_cast<size_t>(j) * stripe_size +
                static_cast<size_t>(seg) * W + lane];
}

// BackTrackBestAlignmentStriped — reads directly from striped SIMD layout.
//
// This is functionally identical to FastPassAligner::BackTrackBestAlignment
// but accesses the M, E, F matrices using the striped index function instead
// of row-major indexing.
SimdGlobalAlignment BackTrackBestAlignmentStriped(
    absl::string_view query, absl::string_view target,
    const int16_t* M_striped, const int16_t* E_striped,
    const int16_t* F_striped, const int16_t* M_row0,
    int S, size_t W, size_t stripe_size,
    int match, int mismatch, int gap_extend) {
  const int n = query.size();
  const int m = target.size();
  const int16_t kNegInf = -30000;

  // Helpers to access matrix values at (i, j).
  auto get_M = [&](int i, int j) -> int {
    if (i == 0) return static_cast<int>(M_row0[j]);
    return static_cast<int>(GetStriped(M_striped, i, j, S, W, stripe_size));
  };
  auto get_E = [&](int i, int j) -> int {
    if (i == 0) return static_cast<int>(kNegInf);
    return static_cast<int>(GetStriped(E_striped, i, j, S, W, stripe_size));
  };
  auto get_F = [&](int i, int j) -> int {
    if (i == 0) return static_cast<int>(kNegInf);
    return static_cast<int>(GetStriped(F_striped, i, j, S, W, stripe_size));
  };

  // Find best score in last row (i=n).
  int best_score = static_cast<int>(kNegInf);
  int best_j = -1;

  for (int j = 0; j <= m; ++j) {
    int score_at_n_j = get_M(n, j);
    if (score_at_n_j >= best_score) {
      best_score = score_at_n_j;
      best_j = j;
    }
  }

  SimdGlobalAlignment alignment;
  alignment.sw_score = best_score;
  alignment.ref_end = best_j - 1;
  alignment.query_end = n - 1;

  int i = n;
  int j = best_j;

  std::vector<char> ops;
  enum Matrix { M_MATRIX, E_MATRIX, F_MATRIX };
  Matrix state = M_MATRIX;

  while (i > 0 || j > 0) {
    if (state == M_MATRIX) {
      if (get_M(i, j) == 0 && (i == 0 || j == 0)) break;
      int score = -60000;  // Very negative sentinel
      if (i > 0 && j > 0) {
        score = (query[i - 1] == target[j - 1]) ? match : mismatch;
      }
      // Highest priority: alignment match/mismatch (diagonal)
      if (i > 0 && j > 0 && get_M(i, j) == get_M(i - 1, j - 1) + score) {
        ops.push_back((query[i - 1] == target[j - 1]) ? '=' : 'X');
        i--;
        j--;
        state = M_MATRIX;
      } else if (i > 0 && get_M(i, j) == get_F(i, j)) {
        state = F_MATRIX;
      } else if (j > 0 && get_M(i, j) == get_E(i, j)) {
        state = E_MATRIX;
      } else {
        break;
      }
    } else if (state == F_MATRIX) {
      ops.push_back('I');
      i--;
      if (i > 0 && get_F(i + 1, j) == get_F(i, j) + gap_extend) {
        state = F_MATRIX;
      } else {
        state = M_MATRIX;
      }
    } else {  // E_MATRIX
      ops.push_back('D');
      j--;
      if (j > 0 && get_E(i, j + 1) == get_E(i, j) + gap_extend) {
        state = E_MATRIX;
      } else {
        state = M_MATRIX;
      }
    }
  }

  std::reverse(ops.begin(), ops.end());
  std::string cigar_str;
  if (!ops.empty()) {
    char last_op = ops[0];
    int count = 1;
    for (size_t k = 1; k < ops.size(); ++k) {
      if (ops[k] == last_op) {
        count++;
      } else {
        absl::StrAppend(&cigar_str, count, std::string(1, last_op));
        last_op = ops[k];
        count = 1;
      }
    }
    absl::StrAppend(&cigar_str, count, std::string(1, last_op));
  }

  alignment.ref_begin = j;
  alignment.query_begin = i;
  if (i > 0) {
    cigar_str = absl::StrCat(i, "S", cigar_str);
  }

  alignment.cigar_string = cigar_str;
  return alignment;
}

}  // namespace

// GlobalAlignSimd — Public entry point
SimdGlobalAlignment GlobalAlignSimd(
    absl::string_view query, absl::string_view target, int edge_range,
    uint8_t match_score, uint8_t mismatch_penalty,
    uint8_t gap_opening_penalty, uint8_t gap_extending_penalty) {
  const int n = query.size();
  const int m = target.size();
  if (n == 0 || m == 0) return SimdGlobalAlignment();

  const int match = static_cast<int>(match_score);
  const int mismatch = -static_cast<int>(mismatch_penalty);
  const int gap_open = -static_cast<int>(gap_opening_penalty);
  const int gap_extend = -static_cast<int>(gap_extending_penalty);
  const int edge_indel_penalty = -static_cast<int>(gap_opening_penalty);

  // Safety check: fall back to scalar path when the DP values may overflow
  // int16. The maximum positive score is query_length * match_score (all bases
  // match), and the most negative score is gap_open + max(n,m) * gap_extend
  // (one long gap). Both must fit within kSimdInt16SafeMax (~30000), which is
  // a conservative limit below int16 max (32767) to leave headroom. With the
  // default match_score of 4, this allows queries up to ~7500 bp.
  constexpr int kSimdInt16SafeMax = 30000;
  const int max_possible_score = n * match;
  const int min_possible_gap = gap_open + std::max(n, m) * gap_extend;
  if (max_possible_score > kSimdInt16SafeMax ||
      min_possible_gap < -kSimdInt16SafeMax) {
    SimdGlobalAlignment overflow_sentinel;
    overflow_sentinel.sw_score = kSimdOverflowSentinel;
    return overflow_sentinel;
  }

  const size_t W = GetSimdWidth();
  const int S = (n + static_cast<int>(W) - 1) / static_cast<int>(W);
  const size_t stripe_size = static_cast<size_t>(S) * W;

  // Allocate full-column striped storage for M, E, F.
  // Each matrix: (m+1) columns × stripe_size values per column.
  const size_t total_striped = static_cast<size_t>(m + 1) * stripe_size;
  auto a_M = hwy::AllocateAligned<int16_t>(total_striped);
  auto a_E = hwy::AllocateAligned<int16_t>(total_striped);
  auto a_F = hwy::AllocateAligned<int16_t>(total_striped);

  // Row 0 storage (separate, since i=0 doesn't map to a query position).
  std::vector<int16_t> M_row0(m + 1);
  std::vector<int16_t> E_row0(m + 1);
  std::vector<int16_t> F_row0(m + 1);

  // Run the SIMD DP.
  HWY_DYNAMIC_DISPATCH(PopulateDpMatrixSimd)
  (query.data(), n, target.data(), m, match, mismatch, gap_open, gap_extend,
   edge_range, edge_indel_penalty, a_M.get(), a_E.get(), a_F.get(),
   M_row0.data(), E_row0.data(), F_row0.data(), stripe_size, S);

  // Run traceback directly on striped layout.
  return BackTrackBestAlignmentStriped(
      query, target, a_M.get(), a_E.get(), a_F.get(), M_row0.data(),
      S, W, stripe_size, match, mismatch, gap_extend);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // HWY_ONCE
