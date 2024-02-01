/*
 * Copyright 2023 Google LLC.
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

#include "deepvariant/alt_aligned_pileup_lib.h"

#include <cstdint>

#include "absl/container/flat_hash_set.h"
#include "third_party/nucleus/protos/cigar.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using CigarUnit = nucleus::genomics::v1::CigarUnit;

namespace {
bool IsOperationRefAdvancing(const CigarUnit& unit) {
  static const absl::flat_hash_set<int>* kRefAdvancingOps =
      new absl::flat_hash_set<int>(
          {nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MATCH,
           nucleus::genomics::v1::CigarUnit::DELETE,
           nucleus::genomics::v1::CigarUnit::SKIP,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MISMATCH});
  return kRefAdvancingOps->contains(unit.operation());
}
bool IsOperationReadAdvancing(const CigarUnit& unit) {
  static const absl::flat_hash_set<int>* kReadAdvancingOps =
      new absl::flat_hash_set<int>(
          {nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MATCH,
           nucleus::genomics::v1::CigarUnit::INSERT,
           nucleus::genomics::v1::CigarUnit::CLIP_SOFT,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MISMATCH});
  return kReadAdvancingOps->contains(unit.operation());
}
}  // namespace

void TrimCigar(const ::google::protobuf::RepeatedPtrField<CigarUnit>& cigar,
               int64_t ref_start, int64_t ref_length,
               ::google::protobuf::RepeatedPtrField<CigarUnit>* new_cigar,
               int64_t& read_start, int64_t& new_read_length) {
  // First consume the ref until the trim is covered.
  int64_t trim_remaining = ref_start;
  // Then consume the ref until the ref_length is covered.
  int64_t ref_to_cover_remaining = ref_length;
  read_start = 0;
  new_read_length = 0;
  for (const auto& cigar_unit : cigar) {
    int64_t c_operation_length = cigar_unit.operation_length();
    // Each operation moves forward in the ref, the read, or both.
    bool advances_ref = IsOperationRefAdvancing(cigar_unit);
    bool advances_read = IsOperationReadAdvancing(cigar_unit);
    int64_t ref_step = advances_ref ? c_operation_length : 0L;
    // First, use up each operation until the trimmed area is covered.
    if (trim_remaining > 0) {
      if (ref_step <= trim_remaining) {
        // Fully apply to the trim.
        trim_remaining -= ref_step;
        read_start += (advances_read ? c_operation_length : 0L);
        continue;
      } else {
        // Partially apply to finish the trim.
        ref_step -= trim_remaining;
        read_start += (advances_read ? trim_remaining : 0L);
        // If trim finishes here, the rest of the ref_step can apply to the
        // next stage and count towards covering the given ref window.
        c_operation_length = ref_step;
        trim_remaining = 0;
      }
    }
    // Once the trim is done, start applying cigar entries to covering the ref
    // window.
    if (trim_remaining == 0) {
      if (ref_step <= ref_to_cover_remaining) {
        // Fully apply to the window.
        // CigarUnit new_cigar_unit;
        auto new_cigar_unit = new_cigar->Add();
        new_cigar_unit->set_operation(cigar_unit.operation());
        new_cigar_unit->set_operation_length(c_operation_length);
        ref_to_cover_remaining -= ref_step;
        new_read_length += (advances_read ? c_operation_length : 0L);
      } else {
        // Partially apply to finish the window.
        c_operation_length = ref_to_cover_remaining;
        auto new_cigar_unit = new_cigar->Add();
        new_cigar_unit->set_operation(cigar_unit.operation());
        new_cigar_unit->set_operation_length(c_operation_length);
        new_read_length += (advances_read ? c_operation_length : 0L);
        ref_to_cover_remaining = 0;
        break;
      }
    }
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
