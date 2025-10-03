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
#include "deepvariant/sampling_util.h"
#include "absl/algorithm/container.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

using nucleus::genomics::v1::Read;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

// Define a tuple type for sorting:
// <hap_index, allele_support_group, position, read_ptr, image_row_ptr>
using ReadPileupTuple =
    std::tuple<int, int, int, const Read*, std::unique_ptr<ImageRow>>;

bool SortImageRows(const ReadPileupTuple& a, const ReadPileupTuple& b) {
  // Primary sort key: haplotype index (std::get<0>(a)).
  if (std::get<0>(a) != std::get<0>(b)) {
    return std::get<0>(a) < std::get<0>(b);
  }

  // Secondary sort key: allele_support_group (std::get<1>(a)).
  // Smaller groups come first.
  if (std::get<1>(a) != std::get<1>(b)) {
    return std::get<1>(a) < std::get<1>(b);
  }

  // Tertiary sort key: alignment position (std::get<2>(a)).
  const Read* read1 = std::get<3>(a);
  const Read* read2 = std::get<3>(b);
  int position1 = std::get<2>(a);
  int position2 = std::get<2>(b);

  if (position1 != position2) {
    return position1 < position2;
  }

  // Tie-breaking: fragment_name, then read_number.
  return std::tuple<std::string, int>(read1->fragment_name(),
                                      read1->read_number()) <
         std::tuple<std::string, int>(read2->fragment_name(),
                                      read2->read_number());
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
    absl::string_view alt_aligned_representation) {
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

std::vector<int> DownsampleReadIndices(
    const std::vector<const ::nucleus::genomics::v1::Read*>& reads,
    int max_reads, std::mt19937_64 gen) {
  // TODO: Use the sampling util function instead.
  std::vector<int> read_indices(reads.size());
  std::iota(read_indices.begin(), read_indices.end(), 0);
  if (reads.size() > max_reads) {
    // Shuffle the indices instead of the reads, so that we won't change the
    // order of the reads list.
    std::shuffle(read_indices.begin(), read_indices.end(), gen);
  }
  return read_indices;
}

AltAlignedPileup GetAltAlignedPileup(std::string alt_aligned_pileup_name) {
  if (alt_aligned_pileup_name == "none") {
    return AltAlignedPileup::kNone;
  } else if (alt_aligned_pileup_name == "base_channels") {
    return AltAlignedPileup::kBaseChannels;
  } else if (alt_aligned_pileup_name == "diff_channels") {
    return AltAlignedPileup::kDiffChannels;
  } else if (alt_aligned_pileup_name == "rows") {
    return AltAlignedPileup::kRows;
  } else if (alt_aligned_pileup_name == "single_row") {
    return AltAlignedPileup::kSingleRow;
  }
  LOG(FATAL) << "Unknown value is specified for alt_aligned_pileup";
  return AltAlignedPileup::kNone;
}

AltAlignedPileup GetSampleAltAlignedPileup(
    AltAlignedPileup global_alt_aligned_pileup,
    std::string sample_alt_aligned_pileup_name) {
  if (!sample_alt_aligned_pileup_name.empty()) {
    return GetAltAlignedPileup(sample_alt_aligned_pileup_name);
  } else {
    return global_alt_aligned_pileup;
  }
}

std::vector<int> GetAltImageRowIndices(
    AltAlignedPileup alt_aligned_pileup,
    absl::Span<const std::string> alt_combination) {
  switch (alt_aligned_pileup) {
    case AltAlignedPileup::kRows:
      return {0, 1};
    case AltAlignedPileup::kSingleRow:
      // If there are 2 alt combinations, return the index of the longer one.
      if (alt_combination.size() == 2 &&
          alt_combination[1].size() > alt_combination[0].size()) {
        return {1};
      }
      return {0};
    default:
      return {};
  }
}

std::vector<int> GetSampleAltImageRowIndices(
    AltAlignedPileup global_alt_aligned_pileup,
    std::string sample_alt_aligned_pileup_name,
    absl::Span<const std::string> alt_combination) {
  AltAlignedPileup sample_alt_aligned_pileup = GetSampleAltAlignedPileup(
      global_alt_aligned_pileup, sample_alt_aligned_pileup_name);
  return GetAltImageRowIndices(sample_alt_aligned_pileup, alt_combination);
}

int CalculatePileupImageHeight(const MakeExamplesOptions& options) {
  int pileup_image_height = 0;
  AltAlignedPileup global_alt_aligned_pileup =
      GetAltAlignedPileup(options.pic_options().alt_aligned_pileup());
  AltAlignedPileup sample_alt_aligned_pileup;
  int num_rows_for_sample;
  for (const auto& sample_options : options.sample_options()) {
    sample_alt_aligned_pileup = GetSampleAltAlignedPileup(
        global_alt_aligned_pileup, sample_options.alt_aligned_pileup());
    // Calculate the number of rows for the sample.
    if (sample_alt_aligned_pileup == AltAlignedPileup::kRows) {
      num_rows_for_sample = 3;
    } else if (sample_alt_aligned_pileup == AltAlignedPileup::kSingleRow) {
      num_rows_for_sample = 2;
    } else {
      num_rows_for_sample = 1;
    }
    pileup_image_height += sample_options.pileup_height() * num_rows_for_sample;
  }
  return pileup_image_height;
}

// Returns a vector of vectors, where each inner vector represents a partition
// of read indices supporting a single allele. The last partition is for reads
// supporting the reference allele.
absl::btree_set<absl::btree_set<int>> GetReadIndicesAllelePartition(
    const DeepVariantCall& dv_call,
    const std::vector<const ::nucleus::genomics::v1::Read*>& reads) {
  // Map read names to their indices in the reads vector. This is used to
  // convert the read names in the DeepVariantCall proto to read indices.
  absl::flat_hash_map<string, int> read_name_to_index_map;
  for (int i = 0; i < reads.size(); ++i) {
    std::string key = (reads[i]->fragment_name() + "/" +
                       std::to_string(reads[i]->read_number()));
    read_name_to_index_map[key] = i;
  }

  // Create a partition element for each allele.
  absl::btree_set<absl::btree_set<int>> read_index_partition_by_allele;

  // Iterate over the allele support map.
  for (const auto& [allele_key, supporting_reads] : dv_call.allele_support()) {
    absl::btree_set<int> read_indices_supporting_allele;
    for (const std::string& read_name : supporting_reads.read_names()) {
      auto it = read_name_to_index_map.find(read_name);
      if (it != read_name_to_index_map.end()) {
        read_indices_supporting_allele.insert(it->second);
        // Remove the read name from the map, so after this section we will only
        // have the reads that don't support any allele.
        read_name_to_index_map.erase(it);
      }
    }
    read_index_partition_by_allele.insert(read_indices_supporting_allele);
  }

  // Ref support info is not always available, so we assume that reads that do
  // not support any allele support the ref.
  absl::btree_set<int> read_indices_supporting_ref;
  for (const auto& [_, index] : read_name_to_index_map) {
    read_indices_supporting_ref.insert(index);
  }

  read_index_partition_by_allele.insert(read_indices_supporting_ref);
  return read_index_partition_by_allele;
}

absl::StatusOr<absl::btree_set<int>> DownsampleReadIndicesWithMinsPerAllele(
    const std::vector<const ::nucleus::genomics::v1::Read*>& reads,
    int max_reads, const DeepVariantCall& dv_call, int min_per_allele,
    std::mt19937_64 gen) {
  absl::btree_set<absl::btree_set<int>> allele_to_read_indices_map =
      GetReadIndicesAllelePartition(dv_call, reads);
  return sampling::SampleWithPartitionMins(allele_to_read_indices_map,
                                           max_reads, min_per_allele, gen);
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
  // TODO Calculate the reference band once and reuse it for each
  // row.
  for (int i = 0; i < options_.reference_band_height(); i++) {
    rows.push_back(EncodeReference(ref_bases));
  }

  // Create a downsampled vector of read indices.
  std::vector<int> sampled_indices(max_reads);
  auto gen = std::mt19937_64(options_.random_seed());
  if (sample_options.use_non_uniform_downsampling()) {
    // Sampling with thresholds may fail if the threshold is too high, so we
    // fall back to uniform sampling if that happens.
    auto status_or_sampled_indices = DownsampleReadIndicesWithMinsPerAllele(
        reads, max_reads, dv_call,
        sample_options.non_uniform_downsampling_threshold(), gen);
    if (!status_or_sampled_indices.ok()) {
      LOG(WARNING) << "Failed to downsample reads with thresholds: "
                   << status_or_sampled_indices.status();
      sampled_indices = DownsampleReadIndices(reads, max_reads, gen);
    } else {
      sampled_indices.assign(status_or_sampled_indices->begin(),
                             status_or_sampled_indices->end());
    }
  } else {
    sampled_indices = DownsampleReadIndices(reads, max_reads, gen);
  }

  // Precompute read-to-allele_group mapping if sorting by alt allele support.
  absl::flat_hash_map<std::string, int> read_name_to_allele_group_map;
  int num_alt_alleles_in_variant = 0;

  if (options_.sort_by_alt_allele_support()) {
    num_alt_alleles_in_variant = dv_call.variant().alternate_bases_size();
    for (int i = 0; i < num_alt_alleles_in_variant; ++i) {
      const std::string& alt = dv_call.variant().alternate_bases(i);
      auto it = dv_call.allele_support().find(alt);
      if (it != dv_call.allele_support().end()) {
        for (const std::string& read_name : it->second.read_names()) {
          read_name_to_allele_group_map[read_name] = i;
        }
      }
    }
  }

  // Each tuple contains:
  // <hap_index, allele_support_group, original_read_alignment_position,
  // read_ptr, image_row_ptr>
  std::vector<ReadPileupTuple> pileup_of_reads;
  for (int index : sampled_indices) {
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

    int hap_idx = GetHapIndex(read);
    int allele_support_group;
    if (options_.sort_by_alt_allele_support()) {
      // Default to a group that sorts after all specific alt alleles.
      allele_support_group = num_alt_alleles_in_variant;
      std::string read_key =
          read.fragment_name() + "/" + std::to_string(read.read_number());
      auto map_it = read_name_to_allele_group_map.find(read_key);
      if (map_it != read_name_to_allele_group_map.end()) {
        allele_support_group = map_it->second;
      }
    } else {
      // If not sorting by allele support, give all reads the same group.
      allele_support_group = 0;
    }

    int64_t read_align_pos =
        (alignment_positions == nullptr || alignment_positions->empty())
            ? read.alignment().position().position()
            : alignment_positions->at(index);

    pileup_of_reads.emplace_back(hap_idx, allele_support_group,
                                 static_cast<int>(read_align_pos), &read,
                                 std::move(image_row));
  }

  absl::c_stable_sort(pileup_of_reads,
                   SortImageRows);
  for (auto& [hap_idx, allele_group, pos, read_ptr, row_ptr] :
       pileup_of_reads) {
    rows.push_back(std::move(row_ptr));
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
    const DeepVariantCall& dv_call, absl::string_view ref_bases,
    const Read& read, int image_start_pos,
    const vector<std::string>& alt_alleles,
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
  img_row.channel_data =
      std::vector<std::vector<unsigned char>>(channel_enums_.size());
  for (int i = 0; i < channel_enums_.size(); ++i) {
    img_row.channel_data[i] = std::vector<unsigned char>(ref_bases.size(), 0);
  }
  bool ok = channel_set.CalculateChannels(
      img_row.channel_data, channel_enums_, read, ref_bases, dv_call,
      alt_alleles, image_start_pos, channels_enum_to_blank);
  // Bail out if we found an issue while calculating channels
  // (a low-quality base at the call site, mapping quality is too low, etc)
  if (!ok) {
    return nullptr;
  }

  return std::make_unique<ImageRow>(img_row);
}

std::unique_ptr<ImageRow> PileupImageEncoderNative::EncodeReference(
    const string& ref_bases) {
  int num_channels = AllChannelsEnum("").size();
  ImageRow img_row(ref_bases.size(), num_channels);
  // Calculate reference rows at the top of each channel image.
  // These are retrieved for each position in the loop below.
  Channels channel_set{options_};
  img_row.channel_data =
      std::vector<std::vector<unsigned char>>(channel_enums_.size());
  for (int i = 0; i < channel_enums_.size(); ++i) {
    img_row.channel_data[i] = std::vector<unsigned char>(ref_bases.size(), 0);
  }
  channel_set.CalculateRefRows(img_row.channel_data, channel_enums_, ref_bases);

  return std::make_unique<ImageRow>(img_row);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
