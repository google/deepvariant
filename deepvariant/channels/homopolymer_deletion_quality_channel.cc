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

#include "deepvariant/channels/homopolymer_deletion_quality_channel.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include "deepvariant/channels/channel.h"

namespace learning {
namespace genomics {
namespace deepvariant {

HomopolymerDeletionQualityChannel::HomopolymerDeletionQualityChannel(
    int width, const PileupImageOptions& options)
    : Channel(width, options) {}

void HomopolymerDeletionQualityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!homopolymer_deletion_quality_vector_.has_value()) {
    homopolymer_deletion_quality_vector_ =
        HomoPolymerInDelQuality(read, true);  // true = deletion
  }
  if (col >= 0 && col < homopolymer_deletion_quality_vector_->size()) {
    data.push_back((*homopolymer_deletion_quality_vector_)[col]);
  } else {
    data.push_back(0);
  }
}

void HomopolymerDeletionQualityChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data.push_back(0);
}

std::vector<int8_t> HomopolymerDeletionQualityChannel::GetTPValues(
    const Read& read) {
  std::vector<int8_t> int_tps(read.aligned_sequence().size());
  if (!read.info().contains("tp")) {
    // Return empty vector if tp tag is not present
    return int_tps;
  }
  const auto& tps = read.info().at("tp").values();

  if (tps.empty()) {
    return int_tps;
  }

  for (int i = 0; i < tps.size(); i++) {
    int_tps[i] = tps[i].int_value();
  }
  return int_tps;
}

std::vector<std::uint8_t>
HomopolymerDeletionQualityChannel::HomoPolymerWeighted(const Read& read) {
  // Generates a vector reflecting the number of repeats observed
  std::vector<std::uint8_t> homopolymer_weighted(read.aligned_sequence().size(),
                                                 1);
  const std::string& seq = read.aligned_sequence();

  if (seq.empty()) {
    return homopolymer_weighted;
  }

  int i = 0;
  while (i < seq.size()) {
    int hmer_length = 1;
    char current_base = seq[i];
    int j = i + 1;

    // Count consecutive identical bases
    while (j < seq.size() && seq[j] == current_base) {
      hmer_length++;
      j++;
    }

    // Fill all positions in this homopolymer with its length
    for (int k = i; k < j; k++) {
      homopolymer_weighted[k] = hmer_length;
    }

    i = j;
  }

  return homopolymer_weighted;
}

std::uint8_t HomopolymerDeletionQualityChannel::BaseQualityColor(
    int base_qual) {
  return static_cast<std::uint8_t>(kMaxPixelValueAsFloat * base_qual /
                                   kMaxQScore);
}

std::vector<std::uint8_t> HomopolymerDeletionQualityChannel::ScaleColorVector(
    std::vector<std::uint8_t>& channel_values, float max_val) {
  std::vector<std::uint8_t> scaled_values(channel_values.size());
  for (int i = 0; i < channel_values.size(); i++) {
    scaled_values[i] = static_cast<std::uint8_t>(kMaxPixelValueAsFloat *
                                                 channel_values[i] / max_val);
  }
  return scaled_values;
}

std::vector<std::uint8_t>
HomopolymerDeletionQualityChannel::HomoPolymerInDelQuality(const Read& read,
                                                           bool is_deletion) {
  // Decode base quality and tp tags to generate a vector indicating the
  // probability of not having an hmer deletion/insertion for each homopolymer
  // in the read. The probabilities are encoded in phred-scores.
  std::vector<std::uint8_t> hmer_directed_qualities(
      read.aligned_sequence().size(), BaseQualityColor(kMaxQScore));

  auto seq = read.aligned_sequence();
  auto hmer_lengths = HomoPolymerWeighted(read);
  auto tps = GetTPValues(read);

  // If tp tag is not present, return default quality values
  if (tps.empty() || tps.size() != seq.size()) {
    return hmer_directed_qualities;
  }

  int i = 0;
  while (i < hmer_lengths.size()) {
    int hmer_length = hmer_lengths[i];
    float hmer_directed_error_prob = 0;

    for (int j = 0; j < hmer_length; j++) {
      if (tps[i + j] == 0) {
        continue;
      }
      bool is_deletion_err = tps[i + j] < 0;
      if (is_deletion_err == is_deletion) {
        std::uint8_t encoded_hmer_qual = read.aligned_quality()[i + j];
        float error_prob = std::pow(10, (encoded_hmer_qual / -10.0));
        hmer_directed_error_prob += error_prob;
      }
    }

    int hmer_directed_quality =
        hmer_directed_error_prob == 0
            ? kMaxQScore
            : static_cast<int>(-10 * std::log10(hmer_directed_error_prob));
    // Clamp to valid range
    if (hmer_directed_quality > kMaxQScore) {
      hmer_directed_quality = kMaxQScore;
    }

    for (int j = 0; j < hmer_length; j++) {
      hmer_directed_qualities[i + j] = BaseQualityColor(hmer_directed_quality);
    }
    i += hmer_length;
  }
  return hmer_directed_qualities;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
