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

#include "deepvariant/channels/read_base_channel.h"

#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

void ReadBaseChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  data[col] = BaseColor(read_base);
}

void ReadBaseChannel::FillRefBase(std::vector<unsigned char>& ref_data, int col,
                                  char ref_base, const std::string& ref_bases) {
  ref_data[col] = BaseColor(ref_base);
}

int ReadBaseChannel::BaseColor(char base) {
  switch (base) {
    case 'A':
      return (options_.base_color_offset_a_and_g() +
              options_.base_color_stride() * 3);
    case 'G':
      return (options_.base_color_offset_a_and_g() +
              options_.base_color_stride() * 2);
    case 'T':
      return (options_.base_color_offset_t_and_c() +
              options_.base_color_stride() * 1);
    case 'C':
      return (options_.base_color_offset_t_and_c() +
              options_.base_color_stride() * 0);
    default:
      return 0;
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
