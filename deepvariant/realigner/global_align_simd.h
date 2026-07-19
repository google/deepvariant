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

// global_align_simd.h — Public API for SIMD-vectorized GlobalAlign
//
// Provides a drop-in replacement for the scalar GlobalAlign that uses
// Google Highway for SIMD-vectorized Farrar striped DP with native
// striped-layout traceback, avoiding the costly layout conversion.

#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_GLOBAL_ALIGN_SIMD_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_GLOBAL_ALIGN_SIMD_H_

#include <cstdint>
#include <string>

#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Result of a SIMD global alignment.
// This mirrors FastPassAligner::GlobalAlignment but is defined independently
// to avoid circular dependencies (fast_pass_aligner depends on this library).
struct SimdGlobalAlignment {
  SimdGlobalAlignment()
      : sw_score(0),
        ref_begin(0),
        ref_end(0),
        query_begin(0),
        query_end(0),
        cigar_string("") {}
  int sw_score;
  int ref_begin;
  int ref_end;
  int query_begin;
  int query_end;
  std::string cigar_string;
};

// Sentinel score indicating SIMD overflow risk — caller should fall back
// to the scalar path.
constexpr int kSimdOverflowSentinel = -2000000000;

// SIMD-vectorized GlobalAlign using Google Highway.
//
// Uses Farrar striped SIMD (int16) for the DP and reads directly from the
// striped layout during traceback, avoiding the O(n*m) conversion that
// negated the speedup in the prototype benchmark.
//
// If query_length * match_score exceeds 30000 (int16 safe limit),
// returns an alignment with sw_score == kSimdOverflowSentinel indicating
// the caller should fall back to the scalar path.
SimdGlobalAlignment GlobalAlignSimd(
    absl::string_view query, absl::string_view target, int edge_range,
    uint8_t match_score, uint8_t mismatch_penalty,
    uint8_t gap_opening_penalty, uint8_t gap_extending_penalty);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_GLOBAL_ALIGN_SIMD_H_
