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

#include "deepvariant/channels/identity_channel.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

IdentityChannel::IdentityChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  identity_color_ = std::nullopt;
}

void IdentityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!identity_color_.has_value()) {
    identity_color_ =
        std::optional<unsigned char>{ScaleColor(Identity(read), kMaxIdentity)};
  }
  data[col] = identity_color_.value();
}

void IdentityChannel::FillRefBase(std::vector<unsigned char>& ref_data, int col,
                                  char ref_base, const std::string& ref_bases) {
  ref_data[col] = static_cast<std::uint8_t>(kMaxPixelValueAsFloat);
}

// Identity: Similar to mapping percent but with a slightly different def.
int IdentityChannel::Identity(const Read& read) {
  int match_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        break;
      case CigarUnit::INSERT:
        break;
      case CigarUnit::DELETE:
        break;
      default:
        break;
    }
  }
  float mapping_percent = (static_cast<float>(match_len) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100;
  return static_cast<int>(mapping_percent);
}

// Scales an input value to pixel range 0-254.
std::uint8_t IdentityChannel::ScaleColor(int value, float max_val) const {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
