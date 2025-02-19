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

#include "deepvariant/channels/avg_base_quality_channel.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/log.h"

namespace learning {
namespace genomics {
namespace deepvariant {

AvgBaseQualityChannel::AvgBaseQualityChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  avg_base_quality_color_ = std::nullopt;
}

void AvgBaseQualityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!avg_base_quality_color_.has_value()) {
    avg_base_quality_color_ = std::optional<unsigned char>{
        ScaleColor(AvgBaseQuality(read), kMaxAvgBaseQuality)};
  }
  data[col] = avg_base_quality_color_.value();
}

void AvgBaseQualityChannel::FillRefBase(std::vector<unsigned char>& ref_data,
                                        int col, char ref_base,
                                        const std::string& ref_bases) {
  ref_data[col] = static_cast<std::uint8_t>(kMaxPixelValueAsFloat);
}

// Scales an input value to pixel range 0-254.
std::uint8_t AvgBaseQualityChannel::ScaleColor(int value, float max_val) const {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}

// Average Base Quality: Averages base quality over length of read.
int AvgBaseQualityChannel::AvgBaseQuality(const Read& read) {
  int base_qual_sum = 0;
  for (const auto& base_qual : read.aligned_quality()) {
    base_qual_sum += base_qual;
    // Base qualities range between 0 and 93
    if (base_qual < 0 || base_qual > 93) {
      LOG(FATAL) << "Encountered base quality outside of bounds (0,93):"
                 << base_qual << ", read=" << read.fragment_name();
    }
  }
  float avg_base_qual = (static_cast<float>(base_qual_sum) /
                         static_cast<float>(read.aligned_quality().size()));
  return static_cast<int>(avg_base_qual);
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
