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

#include "deepvariant/haplotype_tag_channel.h"

#include <cstdint>
#include <string>
#include <vector>

#include "absl/log/log.h"

namespace learning {
namespace genomics {
namespace deepvariant {
void HaplotypeTagChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  const int hp_value =
      GetHPValueForHPChannel(read, options_.hp_tag_for_assembly_polishing());

  read_level_data = std::vector<unsigned char>(1, ScaleColor(hp_value, 2));
}
void HaplotypeTagChannel::FillRefData(const std::string& ref_bases,
                                      std::vector<unsigned char>& ref_data) {
  ref_data = std::vector<unsigned char>(width_, ScaleColor(0, 2));
}

int HaplotypeTagChannel::GetHPValueForHPChannel(
    const Read& read, int hp_tag_for_assembly_polishing) {
  // HP values are added to reads by DeepVariant (direct phasing).
  if (!read.info().contains("HP")) {
    return 0;
  }
  const auto& hp_values = read.info().at("HP").values();
  if (hp_values.empty()) {
    return 0;
  }
  if (hp_values.size() > 1) {
    LOG(WARNING) << "Unexpected: Read contains more than one HP tag. Return 0";
    return 0;
  }
  int hp_value = hp_values[0].int_value();
  // If hp_tag_for_assembly_polishing is set to 2, this is a special case
  // assembly polishing:
  // If we're calling HP=2, when displayed with `--channel_list=haplotype`,
  // we want to swap the color of reads with HP=2 and HP=1.
  if (hp_tag_for_assembly_polishing == 2) {
    if (hp_value == 1) return 2;
    if (hp_value == 2) return 1;
  }

  // Otherwise, keep the default behavior.
  return hp_value;
}

// Scales an input value to pixel range 0-254.
std::uint8_t HaplotypeTagChannel::ScaleColor(int value, float max_val) const {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
