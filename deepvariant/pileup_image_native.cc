/*
 * Copyright 2017 Google LLC.
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

#include "deepvariant/pileup_image_native.h"

#include <math.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "deepvariant/pileup_channel_lib.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

bool SortByAlignment(
    std::tuple<int, int, const Read*, std::unique_ptr<ImageRow>>& a,
    std::tuple<int, int, const Read*, std::unique_ptr<ImageRow>>& b) {
  // Sort reads by position + fragment_name + read_number.
  const Read* read1 = std::get<2>(a);
  const Read* read2 = std::get<2>(b);
  int position1 = std::get<1>(a);
  int position2 = std::get<1>(b);
  if (std::tuple<int, int>(std::get<0>(a), position1) ==
      std::tuple<int, int>(std::get<0>(b), position2)) {
    return std::tuple<std::string, int>(read1->fragment_name(),
                                        read1->read_number()) <
           std::tuple<std::string, int>(read2->fragment_name(),
                                        read2->read_number());
  }
  return std::tuple<int, int>(std::get<0>(a), position1) <
         std::tuple<int, int>(std::get<0>(b), position2);
}

ImageRow::ImageRow(int width, int num_channels)
    : width(width),
      num_channels(num_channels),
      channel_data(num_channels, std::vector<unsigned char>(width, 0)) {}

int ImageRow::Width() const { return width; }

PileupImageEncoderNative::PileupImageEncoderNative(
    const PileupImageOptions& options)
    : options_(options) {
  CHECK((options_.width() % 2 == 1) && options_.width() >= 3)
      << "Width must be odd; found " << options_.width();

  // Pass empty alt_aligned_representation string to AllChannelEnums since
  // alt_aligned_representation is handled at the python level
  channel_enums_ = AllChannelsEnum(/*alt_aligned_representation=*/"");

  // TODO should be CHECK_EQ, but that breaks a test
  CHECK_LE(channel_enums_.size(), options_.num_channels());
}

vector<DeepVariantChannelEnum> PileupImageEncoderNative::AllChannelsEnum(
    const std::string& alt_aligned_representation) {
  std::vector<DeepVariantChannelEnum> channels_list;

  // Fill "default" channels from OptChannel set
  const std::vector<std::string> opt_channels = ToVector(options_.channels());
  for (int j = 0; j < opt_channels.size(); j++) {
    DeepVariantChannelEnum channel =
        Channels::ChannelStrToEnum(opt_channels[j]);
    if (channel != DeepVariantChannelEnum::CH_UNSPECIFIED) {
      channels_list.push_back(channel);
    }
  }

  if (alt_aligned_representation == "diff_channels") {
    channels_list.push_back(
        DeepVariantChannelEnum::CH_DIFF_CHANNELS_ALTERNATE_ALLELE_1);
    channels_list.push_back(
        DeepVariantChannelEnum::CH_DIFF_CHANNELS_ALTERNATE_ALLELE_2);
  } else if (alt_aligned_representation == "base_channels") {
    channels_list.push_back(
        DeepVariantChannelEnum::CH_BASE_CHANNELS_ALTERNATE_ALLELE_1);
    channels_list.push_back(
        DeepVariantChannelEnum::CH_BASE_CHANNELS_ALTERNATE_ALLELE_2);
  }
  return channels_list;
}

std::vector<std::unique_ptr<ImageRow>>
PileupImageEncoderNative::BuildPileupForOneSample(
    const DeepVariantCall& dv_call, const string& ref_bases,
    const std::vector<const ::nucleus::genomics::v1::Read*>& reads,
    int image_start_pos, const vector<std::string>& alt_alleles,
    const SampleOptions& sample_options, const float mean_coverage,
    const std::vector<int64_t>* alignment_positions,
    absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank) {
  // The width of a pileup is defined by the length of ref_bases. ref_bases must
  // have the correct length.
  CHECK(alignment_positions == nullptr || alignment_positions->empty() ||
        alignment_positions->size() == reads.size());
  CHECK_EQ(ref_bases.size(), options_.width());
  int pileup_height = sample_options.pileup_height();
  if (pileup_height == 0) {
    pileup_height = options_.height();
  }
  int max_reads = pileup_height - options_.reference_band_height();

  std::vector<std::unique_ptr<ImageRow>> rows;
  rows.reserve(pileup_height);

  // Reference band at the top of the image.
  for (int i = 0; i < options_.reference_band_height(); i++) {
    rows.push_back(EncodeReference(ref_bases));
  }

  // Create vector of read indices.
  std::vector<int> read_indices(reads.size());
  std::iota(read_indices.begin(), read_indices.end(), 0);
  if (reads.size() > max_reads) {
    // Shuffle the indices instead of the reads, so that we won't change the
    // order of the reads list.
    std::shuffle(read_indices.begin(), read_indices.end(),
                 std::mt19937_64(options_.random_seed()));
  }

  // We add a row for each read in order, down-sampling if the number of
  // reads is greater than the max reads for each sample.
  // Each tuple contains:
  //   <hap_index, original read_alignment_positionm, read, image_row>
  std::vector<std::tuple<int, int, const Read*, std::unique_ptr<ImageRow>>>
      pileup_of_reads;
  for (int index : read_indices) {
    if (pileup_of_reads.size() >= max_reads) {
      break;
    }
    const Read& read = *reads[index];
    std::unique_ptr<ImageRow> image_row =
        EncodeRead(dv_call, ref_bases, read, image_start_pos, alt_alleles,
                  channels_enum_to_blank);
    if (image_row == nullptr) {
      continue;
    }
    if (alignment_positions == nullptr || alignment_positions->empty()) {
      pileup_of_reads.push_back(std::make_tuple(
          GetHapIndex(read), read.alignment().position().position(), &read,
          std::move(image_row)));
    } else {
      pileup_of_reads.push_back(std::make_tuple(GetHapIndex(read),
                                                alignment_positions->at(index),
                                                &read, std::move(image_row)));
    }
  }

  // Sort reads by alignment position.
  std::sort(pileup_of_reads.begin(), pileup_of_reads.end(), SortByAlignment);
  for (auto& [hap_index, pos, read, row] : pileup_of_reads) {
    rows.push_back(std::move(row));
  }

  // Finally, fill in any missing rows to bring our image to pileup_height rows
  // with empty (all black) pixels.
  int empty_rows = pileup_height - rows.size();
  int num_channels = AllChannelsEnum("").size();
  if (empty_rows > 0) {
    for (int i = 0; i < empty_rows; i++) {
      rows.push_back(std::make_unique<ImageRow>(
          ImageRow(ref_bases.size(), num_channels)));
    }
  }

  // Add average coverage information after reads are added and sorted.
  vector<DeepVariantChannelEnum> channels = AllChannelsEnum("");
  auto it = std::find(channel_enums_.begin(), channel_enums_.end(),
                      DeepVariantChannelEnum::CH_MEAN_COVERAGE);
  if (it != channel_enums_.end()) {
    int coverage_channel_index = it - channel_enums_.begin();
    int pileup_width = ref_bases.size();
    for (int i = 0; i < std::min(static_cast<int>(mean_coverage) +
                                     options_.reference_band_height(),
                                 pileup_height);
         i++) {
      if (i < options_.reference_band_height()) {
        // Reference band at the top of the image.
        rows[i]->channel_data[coverage_channel_index].assign(pileup_width,
                                                             kChannelValue255);
      } else {
        // Main part of the image representing mean coverage.
        rows[i]->channel_data[coverage_channel_index].assign(pileup_width,
                                                             kChannelValue200);
      }
    }
  }

  return rows;
}

int PileupImageEncoderNative::GetHapIndex(const Read& read) {
  int default_hap_idx = 0;
  if (!options_.sort_by_haplotypes() || !read.info().contains("HP")) {
    return default_hap_idx;
  }
  nucleus::genomics::v1::ListValue read_info_hp = read.info().at("HP");
  if (read_info_hp.values().empty()) {
    return default_hap_idx;
  }
  const nucleus::genomics::v1::Value& hp_field = read_info_hp.values(0);
  if (hp_field.kind_case() != nucleus::genomics::v1::Value::kIntValue) {
    return default_hap_idx;
  }
  int hp_value = hp_field.int_value();
  int hp_tag_for_assembly_polishing = options_.hp_tag_for_assembly_polishing();
  if (hp_tag_for_assembly_polishing > 0 &&
      hp_value == hp_tag_for_assembly_polishing) {
    // For the target HP tag, set it to -1 so it will be sorted on top of the
    // pileup image.
    return -1;
  } else if (hp_value < 0) {
    // For reads with HP < 0, assume it is not tagged.
    return 0;
  } else {
    return hp_value;
  }
}

std::unique_ptr<ImageRow> PileupImageEncoderNative::EncodeRead(
    const DeepVariantCall& dv_call, const string& ref_bases, const Read& read,
    int image_start_pos, const vector<std::string>& alt_alleles,
    const absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank) {
  int num_channels = AllChannelsEnum("").size();
  ImageRow img_row(ref_bases.size(), num_channels);

  const int mapping_quality = read.alignment().mapping_quality();
  const int min_mapping_quality =
      options_.read_requirements().min_mapping_quality();
  if (mapping_quality < min_mapping_quality) {
    // Bail early if this read's mapping quality is too low.
    return nullptr;
  }

  // Calculate Channels
  Channels channel_set{options_};
  bool ok = channel_set.CalculateChannels(
      channel_enums_, read, ref_bases, dv_call, alt_alleles, image_start_pos,
      channels_enum_to_blank);
  // Bail out if we found an issue while calculating channels
  // (a low-quality base at the call site, mapping quality is too low, etc)
  if (!ok) {
    return nullptr;
  }

  // Fill Channel set.
  img_row.channel_data = std::move(channel_set.data_);
  return std::make_unique<ImageRow>(img_row);
}

std::unique_ptr<ImageRow> PileupImageEncoderNative::EncodeReference(
    const string& ref_bases) {
  int num_channels = AllChannelsEnum("").size();
  ImageRow img_row(ref_bases.size(), num_channels);
  // Calculate reference rows at the top of each channel image.
  // These are retrieved for each position in the loop below.
  Channels channel_set{options_};
  channel_set.CalculateRefRows(channel_enums_, ref_bases);

  img_row.channel_data = std::move(channel_set.ref_data_);

  return std::make_unique<ImageRow>(img_row);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
