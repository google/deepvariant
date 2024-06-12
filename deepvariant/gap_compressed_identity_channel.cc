/*
 * Copyright 2024 Google LLC.
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

#include "deepvariant/gap_compressed_identity_channel.h"

#include <cstdint>
#include <string>
#include <vector>

namespace learning {
namespace genomics {
namespace deepvariant {

void GapCompressedIdentityChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  read_level_data = std::vector<unsigned char>(
      1, ScaleColor(GapCompressedIdentity(read), kMaxIdentity));
}
void GapCompressedIdentityChannel::FillRefData(
    const std::string& ref_bases, std::vector<unsigned char>& ref_data) {
  ref_data = std::vector<unsigned char>(
      width_, static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
}

// Gap Compressed Identity: Ins/Del treated as individual events.
int GapCompressedIdentityChannel::GapCompressedIdentity(const Read& read) {
  // Calculates percentage of the read mapped to the reference
  int match_len = 0;
  int gap_compressed_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        gap_compressed_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        gap_compressed_len += op_len;
        break;
      case CigarUnit::INSERT:
        // Add a single event for insertion.
        gap_compressed_len += 1;
        break;
      case CigarUnit::DELETE:
        // Add a single event for a deletion.
        gap_compressed_len += 1;
        break;
      default:
        break;
    }
  }
  float gap_compressed_identity = static_cast<float>(match_len) /
                                  static_cast<float>(gap_compressed_len) * 100;
  return static_cast<int>(gap_compressed_identity);
}

// Scales an input value to pixel range 0-254.
std::uint8_t GapCompressedIdentityChannel::ScaleColor(int value,
                                                      float max_val) const {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
