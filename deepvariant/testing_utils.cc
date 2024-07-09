/*
 * Copyright 2023 Google LLC.
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

#include "deepvariant/testing_utils.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "google/protobuf/text_format.h"
#include "deepvariant/pileup_image_native.h"
#include "deepvariant/protos/deepvariant.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using ContigInfo = nucleus::genomics::v1::ContigInfo;
using ReferenceSequence = nucleus::genomics::v1::ReferenceSequence;
using Range = nucleus::genomics::v1::Range;

void CreateTestSeq(const std::string& name, const int pos_in_fasta,
                   const int range_start, const int range_end,
                   const std::string& bases, std::vector<ContigInfo>* contigs,
                   std::vector<ReferenceSequence>* seqs) {
  CHECK(pos_in_fasta >= 0 && pos_in_fasta < contigs->size());
  ContigInfo* contig = &contigs->at(pos_in_fasta);
  contig->set_name(name);
  contig->set_pos_in_fasta(pos_in_fasta);
  contig->set_n_bases(range_end - range_start);
  ReferenceSequence* seq = &seqs->at(pos_in_fasta);
  seq->mutable_region()->set_reference_name(name);
  seq->mutable_region()->set_start(range_start);
  seq->mutable_region()->set_end(range_end);
  seq->set_bases(bases);
}

nucleus::genomics::v1::Read MakeRead(
    const absl::string_view chr, const int start, const std::string& bases,
    const std::vector<std::string>& cigar_elements,
    const absl::string_view read_name,
    int hp_tag) {
  nucleus::genomics::v1::Read read = nucleus::MakeRead(
      std::string(chr), start, bases,
      std::vector<std::string>(cigar_elements.begin(), cigar_elements.end()));
  read.set_fragment_name(std::string(read_name));
  if (hp_tag >= 0) {
    auto hp_field = (*read.mutable_info())["HP"].add_values();
    hp_field->set_int_value(hp_tag);
  }
  return read;
}

nucleus::genomics::v1::Variant MakeVariant(
    absl::string_view ref, absl::Span<const absl::string_view> alts,
    int64_t start) {
  nucleus::genomics::v1::Variant variant;
  variant.set_reference_name(std::string(kChr));
  variant.set_start(start);
  variant.set_reference_bases(ref.data(), ref.size());
  for (const auto alt_allele : alts)
    variant.add_alternate_bases(alt_allele.data(), alt_allele.size());

  // End is start + ref length according to Variant.proto spec.
  variant.set_end(variant.start() + ref.length());
  CHECK(google::protobuf::TextFormat::ParseFromString("genotype: -1 genotype: -1",
                                            variant.add_calls()));

  nucleus::genomics::v1::VariantCall* call = variant.mutable_calls(0);
  call->set_call_set_name(std::string(kSampleName));
  return variant;
}

DeepVariantCall MakeDeepVariantCall(
    const nucleus::genomics::v1::Variant& variant) {
  DeepVariantCall dv_call;
  *dv_call.mutable_variant() = variant;
  return dv_call;
}


std::unique_ptr<ImageRow> MakeImageRow(
    const std::vector<std::vector<unsigned char>>& data, int width,
    int num_channels) {
  ImageRow row(width, num_channels);
  int channel_index = 0;
  for (const std::vector<unsigned char>& data_row : data) {
    std::vector<unsigned char> channel_data_row;
    channel_data_row.assign(data_row.begin(), data_row.end());
    row.channel_data[channel_index].assign(data_row.begin(), data_row.end());
    channel_index++;
  }
  return std::make_unique<ImageRow>(row);
}

PileupImageOptions MakeDefaultPileupImageOptions(int width,
                                                 int height,
                                                 int ref_band_height) {
  PileupImageOptions options;
  options.set_reference_band_height(ref_band_height);
  options.set_base_color_offset_a_and_g(40);
  options.set_base_color_offset_t_and_c(30);
  options.set_base_color_stride(70);
  options.set_allele_supporting_read_alpha(1.0);
  options.set_allele_unsupporting_read_alpha(0.6);
  options.set_other_allele_supporting_read_alpha(0.6);
  options.set_reference_matching_read_alpha(0.2);
  options.set_reference_mismatching_read_alpha(1.0);
  options.set_indel_anchoring_base_char("*");
  options.set_reference_alpha(0.4);
  options.set_reference_base_quality(60);
  options.set_positive_strand_color(70);
  options.set_negative_strand_color(240);
  options.set_base_quality_cap(40);
  options.set_mapping_quality_cap(60);
  options.set_height(height);
  options.set_width(width);
  options.set_read_overlap_buffer_bp(5);
  options.set_random_seed(2101079370);
  options.set_min_non_zero_allele_frequency(0.00001);
  return options;
}


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
