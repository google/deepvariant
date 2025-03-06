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

#include "deepvariant/channels/base_6ma_channel.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/io/sam_reader.h"

namespace learning {
namespace genomics {
namespace deepvariant {

Base6mAChannel::Base6mAChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  base_methylation_color_vector_ = std::nullopt;
}

void Base6mAChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!base_methylation_color_vector_.has_value()) {
    std::vector<std::uint8_t> base_modification = GetBaseModification(read);
    base_methylation_color_vector_ = std::optional<std::vector<unsigned char>>{
        ScaleColorVector(base_modification, 255)};
  }
  if (!base_methylation_color_vector_.value().empty()) {
    data[col] = base_methylation_color_vector_.value().at(read_index);
  }
}

void Base6mAChannel::FillRefBase(std::vector<unsigned char>& ref_data,
                                         int col, char ref_base,
                                         const std::string& ref_bases) {
  ref_data[col] = 0;
}

std::vector<std::uint8_t> Base6mAChannel::GetBaseModification(
    const Read& read) {
  auto it = read.base_modifications().find(nucleus::k6mA);
  if (it == read.base_modifications().end()) {
    return {};
  }

  return std::vector<uint8_t>(it->second.begin(), it->second.end());
}

// Scales an input vector to pixel range 0-254
std::vector<std::uint8_t> Base6mAChannel::ScaleColorVector(
    std::vector<std::uint8_t>& channel_values, float max_val) {
  for (int i = 0; i < channel_values.size(); i++) {
    int value = channel_values[i];
    if (static_cast<float>(value) > max_val) {
      value = max_val;
    }
    channel_values[i] = static_cast<int>(kMaxPixelValueAsFloat *
                                         (static_cast<float>(value) / max_val));
  }
  return channel_values;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
