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
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "deepvariant/pileup_channel_lib.h"
#include "deepvariant/protos/deepvariant.pb.h"
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

namespace {

// Get the allele frequency of the alt allele that is carried by a read.
inline float ReadAlleleFrequency(const DeepVariantCall& dv_call,
                                 const Read& read,
                                 const std::vector<std::string>& alt_alleles) {
  return ReadAlleleFrequency_(dv_call, read, alt_alleles);
}

int GetHPValueForHPChannel(const Read& read,
                           int hp_tag_for_assembly_polishing) {
  return GetHPValueForHPChannel_(read, hp_tag_for_assembly_polishing);
}

bool SortByAlignment(std::tuple<int, int, std::unique_ptr<ImageRow>>& a,
                     std::tuple<int, int, std::unique_ptr<ImageRow>>& b) {
  // Sort tuples by making tuples with just the ints.
  return std::tuple<int, int>(std::get<0>(a), std::get<1>(a)) <
         std::tuple<int, int>(std::get<0>(b), std::get<1>(b));
}

}  // namespace

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

// Gets the pixel color (int) for a base.
int PileupImageEncoderNative::BaseColor(char base) const {
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

int PileupImageEncoderNative::BaseColor(const string& base) const {
  CHECK_EQ(base.size(), 1) << "'base' string should be a single character";
  return BaseColor(base[0]);
}

int PileupImageEncoderNative::MatchesRefColor(bool base_matches_ref) const {
  float alpha =
      (base_matches_ref ? options_.reference_matching_read_alpha()
                        : options_.reference_mismatching_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

// Get allele frequency color for a read.
// Convert a frequency value in float to color intensity (int) and normalize.
int PileupImageEncoderNative::AlleleFrequencyColor(
    float allele_frequency) const {
  return AlleleFrequencyColor_(allele_frequency, options_);
}

int PileupImageEncoderNative::SupportsAltColor(int read_supports_alt) const {
  float alpha;
  if (read_supports_alt == 0) {
    alpha = options_.allele_unsupporting_read_alpha();
  } else if (read_supports_alt == 1) {
    alpha = options_.allele_supporting_read_alpha();
  } else {
    CHECK_EQ(read_supports_alt, 2) << "read_supports_alt can only be 0/1/2.";
    alpha = options_.other_allele_supporting_read_alpha();
  }
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

int PileupImageEncoderNative::BaseQualityColor(int base_qual) const {
  float capped =
      static_cast<float>(std::min(options_.base_quality_cap(), base_qual));
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (capped / options_.base_quality_cap()));
}

int PileupImageEncoderNative::MappingQualityColor(int mapping_qual) const {
  float capped = static_cast<float>(
      std::min(options_.mapping_quality_cap(), mapping_qual));
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (capped / options_.mapping_quality_cap()));
}

int PileupImageEncoderNative::StrandColor(bool on_positive_strand) const {
  return (on_positive_strand ? options_.positive_strand_color()
                             : options_.negative_strand_color());
}

vector<DeepVariantChannelEnum> PileupImageEncoderNative::AllChannelsEnum(
    const std::string& alt_aligned_representation) {
  std::vector<DeepVariantChannelEnum> channels_list;
  // The list here corresponds to the order in EncodeRead.
  channels_list.push_back(DeepVariantChannelEnum::CH_READ_BASE);
  channels_list.push_back(DeepVariantChannelEnum::CH_BASE_QUALITY);
  channels_list.push_back(DeepVariantChannelEnum::CH_MAPPING_QUALITY);
  channels_list.push_back(DeepVariantChannelEnum::CH_STRAND);
  channels_list.push_back(DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT);
  channels_list.push_back(DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF);
  if (options_.use_allele_frequency()) {
    channels_list.push_back(DeepVariantChannelEnum::CH_ALLELE_FREQUENCY);
  }
  if (options_.add_hp_channel()) {
    channels_list.push_back(DeepVariantChannelEnum::CH_HAPLOTYPE_TAG);
  }
  // Fill OptChannel set.
  const std::vector<std::string> opt_channels = ToVector(options_.channels());
  for (int j = 0; j < opt_channels.size(); j++) {
    channels_list.push_back(ChannelStrToEnum(opt_channels[j]));
  }
  // Then, in pileup_image.py, _represent_alt_aligned_pileups can potentially
  // add two more channels.
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
    const std::vector<const ::nucleus::genomics::v1::Read *>& reads,
    int image_start_pos,
    const vector<std::string>& alt_alleles,
    const SampleOptions& sample_options) {
  // The width of a pileup is defined by the length of ref_bases. ref_bases must
  // have the correct length.
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
  std::vector<std::tuple<int, int, std::unique_ptr<ImageRow>>> pileup_of_reads;
  for (int index : read_indices) {
    if (pileup_of_reads.size() >= max_reads) {
      break;
    }
    const Read& read = *reads[index];
    std::unique_ptr<ImageRow> image_row =
        EncodeRead(dv_call, ref_bases, read, image_start_pos, alt_alleles);
    if (image_row == nullptr) {
      continue;
    }
    pileup_of_reads.push_back(std::make_tuple(
        GetHapIndex(read), read.alignment().position().position(),
        std::move(image_row)));
  }

  // Sort reads by alignment position.
  std::sort(pileup_of_reads.begin(), pileup_of_reads.end(), SortByAlignment);
  for (auto& [hap_index, position, row] : pileup_of_reads) {
    rows.push_back(std::move(row));
  }

  // Finally, fill in any missing rows to bring our image to pileup_height rows
  // with empty (all black) pixels.
  int empty_rows = pileup_height - rows.size();
  if (empty_rows > 0) {
    for (int i = 0; i < empty_rows; i++) {
      rows.push_back(std::make_unique<ImageRow>(
          ImageRow(ref_bases.size(), options_.num_channels())));
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
  nucleus::genomics::v1::Value hp_field = read_info_hp.values(0);
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
    int image_start_pos, const vector<std::string>& alt_alleles) {
  ImageRow img_row(ref_bases.size(), options_.num_channels());

  const int mapping_quality = read.alignment().mapping_quality();
  const int min_mapping_quality =
      options_.read_requirements().min_mapping_quality();
  // Bail early if this read's mapping quality is too low.
  if (mapping_quality < min_mapping_quality) {
    return nullptr;
  }
  // Calculate Channels
  Channels channel_set{options_};
  bool ok = channel_set.CalculateChannels(
      channel_enums_, read, ref_bases, dv_call, alt_alleles, image_start_pos);
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
  ImageRow img_row(ref_bases.size(), options_.num_channels());

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
