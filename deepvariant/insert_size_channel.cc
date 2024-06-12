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

#include "deepvariant/insert_size_channel.h"

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace learning {
namespace genomics {
namespace deepvariant {

void InsertSizeChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  read_level_data = ReadInsertSize(read);
}
void InsertSizeChannel::FillRefData(const std::string& ref_bases,
                                    std::vector<unsigned char>& ref_data) {
  ref_data = std::vector<unsigned char>(
      width_, static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
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
