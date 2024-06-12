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

#include "deepvariant/gc_content_channel.h"

#include <cstdint>
#include <string>
#include <vector>

namespace learning {
namespace genomics {
namespace deepvariant {
void GcContentChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  read_level_data =
      std::vector<unsigned char>(1, ScaleColor(GcContent(read), kMaxGcContent));
}
void GcContentChannel::FillRefData(const std::string& ref_bases,
                                   std::vector<unsigned char>& ref_data) {
  Read refRead;
  refRead.set_aligned_sequence(ref_bases);
  ref_data = std::vector<unsigned char>(
      width_, ScaleColor(GcContent(refRead), kMaxGcContent));
}

int GcContentChannel::GcContent(const Read& read) {
  int gc_count{};

  for (const auto& base : read.aligned_sequence()) {
    if (base == 'G' || base == 'C') {
      gc_count += 1;
    }
  }

  return static_cast<int>((static_cast<float>(gc_count) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100);
}

// Scales an input value to pixel range 0-254.
std::uint8_t GcContentChannel::ScaleColor(int value, float max_val) const {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
