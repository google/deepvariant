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

#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/channels/homopolymer_indel_quality_channel.h"
#include "deepvariant/protos/deepvariant.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

HomopolymerDeletionQualityChannel::HomopolymerDeletionQualityChannel(
    int width, const PileupImageOptions& options)
    : HomopolymerInDelQualityChannel(width, options) {}

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

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
