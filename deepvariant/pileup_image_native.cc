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
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "deepvariant/pileup_channel_lib.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "absl/log/log.h"

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

using tensorflow::uint8;

namespace {

// Get the allele frequency of the alt allele that is carried by a read.
inline float ReadAlleleFrequency(const DeepVariantCall& dv_call,
                                 const Read& read,
                                 const std::vector<std::string>& alt_alleles) {
  string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
  for (const string& alt_allele : dv_call.variant().alternate_bases()) {
    const auto& allele_support = dv_call.allele_support();
    auto it_read = allele_support.find(alt_allele);

    if (it_read != allele_support.end()) {
      const auto& supp_read_names = it_read->second.read_names();
      for (const string& read_name : supp_read_names) {
        const bool alt_in_alt_alleles =
            std::find(alt_alleles.begin(), alt_alleles.end(), alt_allele) !=
            alt_alleles.end();
        // If the read supports an alt we are currently considering, return the
        // associated allele frequency.
        if (read_name == key && alt_in_alt_alleles) {
          auto it = dv_call.allele_frequency().find(alt_allele);
          if (it != dv_call.allele_frequency().end())
            return it->second;
          else
            return 0;
        }
      }
    }
  }
  // If cannot find the matching variant, set the frequency to 0.
  return 0;
}

int GetHPValueForHPChannel(const Read& read,
                           int hp_tag_for_assembly_polishing) {
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
  // See the description of --add_hp_channel flag in make_examples.py: Currently
  // we only support value of 1, 2, or 0.
  if (hp_value != 0 && hp_value != 1 && hp_value != 2) {
    LOG(FATAL)
        << "This function is currently used when --add_hp_channel is set. "
        << "HP value has to be either 1, 2, or 0. Found a read with HP="
        << hp_value << ", read=" << read.DebugString();
  }
  // If hp_tag_for_assembly_polishing is set to 2, this is a special case
  // assembly polishing:
  // If we're calling HP=2, when displayed with --add_hp_channel, we want to
  // swap the color of reads with HP=2 and HP=1.
  if (hp_tag_for_assembly_polishing == 2) {
    if (hp_value == 1) return 2;
    if (hp_value == 2) return 1;
  }

  // Otherwise, keep the default behavior.
  return hp_value;
}

}  // namespace

ImageRow::ImageRow(int width, int num_channels, bool use_allele_frequency,
                   bool add_hp_channel, const std::vector<string>& channels)
    : base(width, 0),
      base_quality(width, 0),
      mapping_quality(width, 0),
      on_positive_strand(width, 0),
      supports_alt(width, 0),
      matches_ref(width, 0),
      sequencing_type(width, 0),
      allele_frequency(width, 0),
      hp_value(width, 0),
      num_channels(num_channels),
      use_allele_frequency(use_allele_frequency),
      add_hp_channel(add_hp_channel),
      channels(channels) {}

int ImageRow::Width() const {
  CHECK(
      base.size() == base_quality.size() &&
      base.size() == mapping_quality.size() &&
      base.size() == on_positive_strand.size() &&
      base.size() == supports_alt.size() && base.size() == matches_ref.size() &&
      base.size() == sequencing_type.size() &&
      base.size() == allele_frequency.size() && base.size() == hp_value.size());
  return base.size();
}

PileupImageEncoderNative::PileupImageEncoderNative(
    const PileupImageOptions& options)
    : options_(options) {
  CHECK((options_.width() % 2 == 1) && options_.width() >= 3)
      << "Width must be odd; found " << options_.width();
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
  if (allele_frequency <= options_.min_non_zero_allele_frequency()) {
    return 0;
  } else {
    float log10_af = log10(allele_frequency);
    float log10_min = log10(options_.min_non_zero_allele_frequency());
    return ((log10_min - log10_af) / log10_min) *
           static_cast<int>(kMaxPixelValueAsFloat);
  }
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


std::unique_ptr<ImageRow> PileupImageEncoderNative::EncodeRead(
    const DeepVariantCall& dv_call, const string& ref_bases, const Read& read,
    int image_start_pos, const vector<std::string>& alt_alleles) {
  ImageRow img_row(ref_bases.size(), options_.num_channels(),
                   options_.use_allele_frequency(), options_.add_hp_channel(),
                   ToVector(options_.channels()));

  // Calculate base channels.
  const int supports_alt = ReadSupportsAlt(dv_call, read, alt_alleles);
  const int mapping_quality = read.alignment().mapping_quality();
  const int min_mapping_quality =
      options_.read_requirements().min_mapping_quality();
  // Bail early if this read's mapping quality is too low.
  if (mapping_quality < min_mapping_quality) {
    return nullptr;
  }
  const bool is_forward_strand = !read.alignment().position().reverse_strand();
  const uint8 alt_color = SupportsAltColor(supports_alt);
  const uint8 mapping_color = MappingQualityColor(mapping_quality);
  const uint8 strand_color = StrandColor(is_forward_strand);
  const int min_base_quality = options_.read_requirements().min_base_quality();

  // Calculate AUX channels.
  const float allele_frequency =
      (options_.use_allele_frequency())
          ? ReadAlleleFrequency(dv_call, read, alt_alleles)
          : 0;
  const uint8 allele_frequency_color = AlleleFrequencyColor(allele_frequency);
  const int hp_value = (options_.add_hp_channel())
                           ? GetHPValueForHPChannel(
                                 read, options_.hp_tag_for_assembly_polishing())
                           : 0;

  // Calculate OptChannels.
  OptChannels channel_set{options_};
  bool ok = channel_set.CalculateChannels(img_row.channels, read,
                                          ref_bases, dv_call, alt_alleles,
                                          image_start_pos);
  // Bail out if we found an issue while calculating channels
  // (a low-quality base at the call site, mapping quality is too low, etc)
  if (!ok) {
    return nullptr;
  }

  // Fill OptChannel set.
  img_row.channel_data.resize(img_row.channels.size(),
                            std::vector<unsigned char>(ref_bases.size(), 0));
  for (int j = 0; j < img_row.channels.size(); j++) {
    const std::string channel = img_row.channels[j];
    img_row.channel_data[j] = channel_set.data_[channel];
  }

  // Handler for each component of the CIGAR string, as subdivided
  // according the rules below.
  // Side effect: draws in img_row
  // Return value: true on normal exit; false if we determine that we
  // have a low quality base at the call position (in which case we
  // should return null) from EncodeRead.
  std::function<bool(int, int, const CigarUnit::Operation&)>
      action_per_cigar_unit =
          [&](int ref_i, int read_i, const CigarUnit::Operation& cigar_op) {
            char read_base = 0;
            if (cigar_op == CigarUnit::INSERT) {
              // TODO: fix this to be a char in the proto
              read_base = options_.indel_anchoring_base_char()[0];
            } else if (cigar_op == CigarUnit::DELETE) {
              ref_i -= 1;  // Adjust anchor base on reference
              read_base = options_.indel_anchoring_base_char()[0];
            } else if (cigar_op == CigarUnit::ALIGNMENT_MATCH ||
                       cigar_op == CigarUnit::SEQUENCE_MATCH ||
                       cigar_op == CigarUnit::SEQUENCE_MISMATCH) {
              read_base = read.aligned_sequence()[read_i];
            }

            size_t col = ref_i - image_start_pos;
            if (read_base && 0 <= col && col < ref_bases.size()) {
              int base_quality = read.aligned_quality(read_i);
              if (ref_i == dv_call.variant().start() &&
                  base_quality < min_base_quality) {
                return false;
              }
              bool matches_ref = (read_base == ref_bases[col]);

              // Fill Base channel set.
              img_row.base[col] = BaseColor(read_base);
              img_row.base_quality[col] = BaseQualityColor(base_quality);
              img_row.mapping_quality[col] = mapping_color;
              img_row.on_positive_strand[col] = strand_color;
              img_row.supports_alt[col] = alt_color;
              img_row.matches_ref[col] = MatchesRefColor(matches_ref);

              // Fill AUX channel set.
              if (img_row.use_allele_frequency) {
                img_row.allele_frequency[col] = allele_frequency_color;
              }
              if (img_row.add_hp_channel) {
                img_row.hp_value[col] = ScaleColor(hp_value, 2);
              }
            }
            return true;
          };

  // In the following, we iterate over alignment information for each
  // base of read, invoking action_per_cigar_unit for every segment of
  // the alignment.
  //
  // The handling of each cigar element type is given below, assuming
  // it has length n.
  //
  // ALIGNMENT_MATCH, SEQUENCE_MATCH, SEQUENCE_MISMATCH:
  //   Provide a segment ref_i, read_i for each of the n bases in the
  //   operator, where ref_i is the position on the genome where this
  //   base aligns.
  //
  // INSERT, CLIP_SOFT:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the preceding base of the insertion; i.e., the anchor
  //   base. Beware that ref_i could be -1 if the insertion is aligned
  //   to the first base of a contig.  read_i points to the first base
  //   of the insertion. So if our cigar is 1M2I1M for a read starting
  //   at S, we'd see first (S, 0, '1M'), followed by one (S, 1,
  //   '2I'), and then (S + 1, 3, '1M').
  //
  // DELETE, SKIP:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the first base of the deletion, just like in an
  //   ALIGNMENT_MATCH. read_i points to the previous base in the
  //   read, as there's no actual read sequence associated with a
  //   deletion. Beware that read_i could be -1 if the deletion is the
  //   first cigar of the read. So if our cigar is 1M2D1M for a read
  //   starting at S, we'd see first (S, 0, '1M'), followed by one (S
  //   + 1, 0, '2D'), and then (S + 3, 1, '1M').
  //
  // CLIP_HARD, PAD:
  //   These operators are ignored by as they don't impact the
  //   alignment of the read w.r.t. the reference.
  //
  // Any other CIGAR op:
  //   Fatal error, at present; later we should fail with a status encoding.
  int ref_i = read.alignment().position().position();
  int read_i = 0;
  ok = true;

  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();

    switch (op) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
        // Alignment op.
        for (int i = 0; i < op_len; i++) {
          ok = ok && action_per_cigar_unit(ref_i, read_i, op);
          ref_i++;
          read_i++;
        }
        break;
      case CigarUnit::INSERT:
      case CigarUnit::CLIP_SOFT:
        // Insert op.
        ok = action_per_cigar_unit(ref_i - 1, read_i, op);
        read_i += op_len;
        break;
      case CigarUnit::DELETE:
      case CigarUnit::SKIP:
        // Delete op.
        ok = action_per_cigar_unit(ref_i, read_i - 1, op);
        ref_i += op_len;
        break;
      case CigarUnit::CLIP_HARD:
      case CigarUnit::PAD:
        // Ignored ops.  Do nothing.
        break;
      default:
        LOG(FATAL) << "Unrecognized CIGAR op";
    }

    // Bail out if we found this read had a low-quality base at the
    // call site.
    if (!ok) {
      return nullptr;
    }
  }

  return std::make_unique<ImageRow>(img_row);
}

std::unique_ptr<ImageRow> PileupImageEncoderNative::EncodeReference(
    const string& ref_bases) {
  int ref_qual = options_.reference_base_quality();
  uint8 base_quality_color = BaseQualityColor(ref_qual);
  uint8 mapping_quality_color = MappingQualityColor(ref_qual);
  // We use "+" strand color for the reference.
  uint8 strand_color = StrandColor(true);
  uint8 alt_color = SupportsAltColor(0);
  uint8 ref_color = MatchesRefColor(true);
  uint8 allele_frequency_color = AlleleFrequencyColor(0);

  ImageRow img_row(ref_bases.size(), options_.num_channels(),
                   options_.use_allele_frequency(), options_.add_hp_channel(),
                   ToVector(options_.channels()));

  // Initialize dynamic channels
  img_row.channel_data.resize(img_row.channels.size(),
                              std::vector<unsigned char>(ref_bases.size(), 0));

  // Calculate reference rows at the top of each channel image.
  // These are retrieved for each position in the loop below.
  OptChannels channel_set{options_};
  channel_set.CalculateRefRows(img_row.channels, ref_bases);

  for (size_t col = 0; col < ref_bases.size(); ++col) {
    img_row.base[col] = BaseColor(ref_bases[col]);
    img_row.base_quality[col] = base_quality_color;
    img_row.mapping_quality[col] = mapping_quality_color;
    img_row.on_positive_strand[col] = strand_color;
    img_row.supports_alt[col] = alt_color;
    img_row.matches_ref[col] = ref_color;
    if (img_row.use_allele_frequency) {
      img_row.allele_frequency[col] = allele_frequency_color;
    }
    if (img_row.add_hp_channel) {
      img_row.hp_value[col] = ScaleColor(0, 2);
    }

    // Optional channels
    for (int j = 0; j < img_row.channels.size(); j++) {
      img_row.channel_data[j][col] =
          channel_set.GetRefRows(img_row.channels[j], col);
    }
  }

  return std::make_unique<ImageRow>(img_row);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
