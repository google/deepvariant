/*
 * Copyright 2017 Google LLC.
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
#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_WINDOW_SELECTOR_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_WINDOW_SELECTOR_H_

#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/realigner.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Returns a vector of candidate counts for the realigner/window_selector.
//
// This algorithm assumes that the provided AlleleCounter object has had the
// reads we want to use to compute our candidate positions already added to it.
//
// We walk over the AlleleCounts in the AlleleCounter and increment our position
// counts for each non-REFERENCE allele we see as follows:
//
//  - SUBSTITITUION: at the substitution position only.
//  - DELETE: at positions within [position, position + length)
//  - INSERT and CLIP_SOFT: at positions within
//      [position - cigar_len, position + cigar_len)
//
// Returns:
//   A vector of counts. The value at offset [i] is the sum of candidate
//   positions across all alleles in AlleleCount that affected position i.
//   The i-th position in this vector corresponding to the i-th position in
//   allele_count->Interval(). So result[0] is the count for
//   allele_count->Interval().start().
std::vector<int> VariantReadsWindowSelectorCandidates(
    const AlleleCounter& allele_counter);

std::vector<float> AlleleCountLinearWindowSelectorCandidates(
    const AlleleCounter& allele_counter,
    const WindowSelectorModel::AlleleCountLinearModel& config);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_WINDOW_SELECTOR_H_
