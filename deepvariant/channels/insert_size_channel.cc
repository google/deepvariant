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

#include "deepvariant/channels/insert_size_channel.h"

#include <cstdint>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

InsertSizeChannel::InsertSizeChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  insert_size_color_vector_ = std::nullopt;
}

void InsertSizeChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!insert_size_color_vector_.has_value()) {
    insert_size_color_vector_ =
        std::optional<std::vector<unsigned char>>{ReadInsertSize(read)};
  }
  data[col] = insert_size_color_vector_.value().at(read_index);
}

void InsertSizeChannel::FillRefBase(std::vector<unsigned char>& ref_data,
                                    int col, char ref_base,
                                    const std::string& ref_bases) {
  ref_data[col] = static_cast<std::uint8_t>(kMaxPixelValueAsFloat);
}

std::vector<std::uint8_t> InsertSizeChannel::ReadInsertSize(const Read& read) {
  // Generates a vector reflecting the fragment length of the read
  std::vector<std::uint8_t> reads_with_insert_size(
      read.aligned_sequence().size(), normalizeFragmentLength(read));
  return reads_with_insert_size;
}

// normalizes a Read's `fragment_length` to a pixel value
int InsertSizeChannel::normalizeFragmentLength(const Read& read) {
  int fragment_length = std::abs(read.fragment_length());
  if (static_cast<float>(fragment_length) > kMaxFragmentLength) {
    fragment_length = static_cast<int>(kMaxFragmentLength);
  }
  return static_cast<int>(
      kMaxPixelValueAsFloat *
      (static_cast<float>(fragment_length) / kMaxFragmentLength));
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
