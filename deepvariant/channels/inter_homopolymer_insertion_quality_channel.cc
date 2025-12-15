/*
 * Copyright 2025 Google LLC.
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

#include "deepvariant/channels/inter_homopolymer_insertion_quality_channel.h"

#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/channels/channel_utils.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

InterHomopolymerInsertionQualityChannel::
    InterHomopolymerInsertionQualityChannel(int width,
                                            const PileupImageOptions& options)
    : Channel(width, options) {}

void InterHomopolymerInsertionQualityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!inter_homopolymer_insertion_quality_vector_.has_value()) {
    inter_homopolymer_insertion_quality_vector_ = GetT0QualityValues(read);
  }
  if (read_index >= 0 &&
      read_index < inter_homopolymer_insertion_quality_vector_->size()) {
    data[col] = (*inter_homopolymer_insertion_quality_vector_)[read_index];
  } else {
    data[col] = 0;
  }
}

void InterHomopolymerInsertionQualityChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data.push_back(0);
}

std::vector<std::uint8_t> InterHomopolymerInsertionQualityChannel::GetT0Values(
    const Read& read) {
  /*
  Load T0 values for reads.
  T0 values are phred qualities which encode the probability of having a
  non-homopolymeric insertion between every two homopolymers in the read.
  For example:
  For sequence AAAATTTT t0:Z:5555IIII means that there is a probability
  corresponding to Q-score 20 (encoded as '5' in ASCII - 33 = 20) for inserted
  bases (not A*T*) between the homopolymer AAAA and TTTT.

  Return - t0 values as Q-scores
  */
  std::vector<std::uint8_t> int_t0s(read.aligned_sequence().size(), 0);
  if (!read.info().contains("t0")) {
    // Return zero vector if t0 tag is not present.
    return int_t0s;
  }
  const auto& t0_values = read.info().at("t0").values();
  if (t0_values.empty()) {
    return int_t0s;
  }
  if (t0_values.size() > 1) {
    LOG(WARNING) << "Read contains more than one t0 tag.";
  }

  const absl::string_view t0 = t0_values[0].string_value();

  if (t0.empty()) {
    return int_t0s;
  }

  for (int i = 0; i < t0.size() && i < int_t0s.size(); i++) {
    // Convert from ASCII-encoded phred score (char - 33)
    int_t0s[i] = static_cast<std::uint8_t>(t0[i] - 33);
  }
  return int_t0s;
}

std::vector<std::uint8_t>
InterHomopolymerInsertionQualityChannel::GetT0QualityValues(const Read& read) {
  // Get T0 values and convert them to color values
  auto t0_values = GetT0Values(read);
  std::vector<std::uint8_t> t0_quality_colors(t0_values.size());

  for (int i = 0; i < t0_values.size(); i++) {
    t0_quality_colors[i] = channels::internal::BaseQualityColor(t0_values[i]);
  }

  return t0_quality_colors;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
