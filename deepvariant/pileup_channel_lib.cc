

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
#include "deepvariant/pileup_channel_lib.h"

#include <math.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

bool Channels::CalculateChannels(
    const std::vector<DeepVariantChannelEnum>& channel_enums, const Read& read,
    absl::string_view ref_bases, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles, int image_start_pos,
    const absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank) {
  int maxEnumValue = getMaxEnumValue(channel_enums);
  channel_enum_to_index_ = std::vector<int>(maxEnumValue + 1);
  int currIndex = 0;
  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    channel_enum_to_index_[channel_enum] = currIndex;
    currIndex++;
  }

  /*--------------------------------------
  Calculate read-level channels
  ---------------------------------------*/
  data_ = std::vector<std::vector<unsigned char>>(channel_enums.size());
  read_level_data_ =
      std::vector<std::vector<unsigned char>>(channel_enums.size());

  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    int index = channel_enum_to_index_[channel_enum];
    data_[index] = std::vector<unsigned char>(ref_bases.size(), 0);
    if (!isBaseLevelChannel(channel_enum) &&
        !channels_enum_to_blank.contains(channel_enum)) {
      bool ok =
          CalculateReadLevelData(channel_enum, read, dv_call, alt_alleles);
      if (!ok) return false;
    }
  }
  /*--------------------------------------
  Calculate base-level channels
  ---------------------------------------*/

  // Handler for each component of the CIGAR string, as subdivided
  // according the rules below.
  // Side effect: draws in img_row
  // Return value: true on normal exit; false if we determine that we
  // have a low quality base at the call position (in which case we
  // should return null) from EncodeRead.
  std::function<bool(int, int, const CigarUnit::Operation&)>
      action_per_cigar_unit = [&, ref_bases](
                                  int ref_i, int read_i,
                                  const CigarUnit::Operation& cigar_op) {
        char read_base = 0;
        if (cigar_op == CigarUnit::INSERT) {
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
          // Bail out if we found this read had a low-quality base at the
          // call site.
          if (ref_i == dv_call.variant().start() &&
              base_quality < options_.read_requirements().min_base_quality()) {
            return false;
          }

          // Calculate base level values for channels
          for (int i = 0; i < channel_enums.size(); ++i) {
            DeepVariantChannelEnum channel_enum = channel_enums[i];
            int index = channel_enum_to_index_[channel_enum];
            if (isBaseLevelChannel(channel_enum) &&
                !channels_enum_to_blank.contains(channel_enum)) {
              if (channel_enum == DeepVariantChannelEnum::CH_READ_BASE) {
                data_[index][col] = BaseColor(read_base, options_);
              } else if (channel_enum ==
                         DeepVariantChannelEnum::CH_BASE_QUALITY) {
                data_[index][col] =
                    ScaleColor(base_quality, options_.base_quality_cap());
              } else if (channel_enum ==
                         DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF) {
                bool matches_ref = (read_base == ref_bases[col]);
                data_[index][col] = MatchesRefColor(matches_ref, options_);
              }
            }
          }

          // Fill in base level value for read level channels from
          // previously calculated read level values
          for (int i = 0; i < channel_enums.size(); ++i) {
            DeepVariantChannelEnum channel_enum = channel_enums[i];
            int index = channel_enum_to_index_[channel_enum];
            if (!isBaseLevelChannel(channel_enum) &&
                !channels_enum_to_blank.contains(channel_enum)) {
              if (read_level_data_[index].size() == 1) {
                data_[index][col] = read_level_data_[index][0];
              } else {
                data_[index][col] = read_level_data_[index][read_i];
              }
            }
          }
        }
        return true;
      };
  return CalculateBaseLevelData(read, action_per_cigar_unit);
}

// Calculate values for channels that only depend on information at the
// granularity of an entire read.
bool Channels::CalculateReadLevelData(
    const DeepVariantChannelEnum channel_enum, const Read& read,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  int index = channel_enum_to_index_[channel_enum];

  if (channel_enum == DeepVariantChannelEnum::CH_HAPLOTYPE_TAG) {
    const int hp_value =
        GetHPValueForHPChannel(read, options_.hp_tag_for_assembly_polishing());
    read_level_data_[index].assign({ScaleColor(hp_value, 2)});
  } else if (channel_enum == DeepVariantChannelEnum::CH_ALLELE_FREQUENCY) {
    const float allele_frequency =
        ReadAlleleFrequency(dv_call, read, alt_alleles);
    read_level_data_[index].assign({static_cast<std::uint8_t>(
        AlleleFrequencyColor(allele_frequency, options_))});
  } else if (channel_enum == DeepVariantChannelEnum::CH_MAPPING_QUALITY) {
    const int mapping_quality = read.alignment().mapping_quality();
    const int min_mapping_quality =
        options_.read_requirements().min_mapping_quality();
    // Bail early if this read's mapping quality is too low.
    if (mapping_quality < min_mapping_quality) {
      return false;
    }
    read_level_data_[index].assign(
        {ScaleColor(mapping_quality, options_.mapping_quality_cap())});
  } else if (channel_enum == DeepVariantChannelEnum::CH_STRAND) {
    const bool is_forward_strand =
        !read.alignment().position().reverse_strand();
    read_level_data_[index].assign(
        {static_cast<std::uint8_t>(StrandColor(is_forward_strand, options_))});
  } else if (channel_enum == DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT) {
    int supports_alt = ReadSupportsAlt(dv_call, read, alt_alleles);
    read_level_data_[index].assign(
        {static_cast<std::uint8_t>(SupportsAltColor(supports_alt, options_))});
  } else if (channel_enum == DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT) {
    read_level_data_[index].assign(
        {ScaleColor(ReadMappingPercent(read), MaxMappingPercent)});
  } else if (channel_enum == DeepVariantChannelEnum::CH_AVG_BASE_QUALITY) {
    read_level_data_[index].assign(
        {ScaleColor(AvgBaseQuality(read), MaxAvgBaseQuality)});
  } else if (channel_enum == DeepVariantChannelEnum::CH_IDENTITY) {
    read_level_data_[index].assign({ScaleColor(Identity(read), MaxIdentity)});
  } else if (channel_enum ==
             DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY) {
    read_level_data_[index].assign(
        {ScaleColor(GapCompressedIdentity(read), MaxIdentity)});
  } else if (channel_enum == DeepVariantChannelEnum::CH_GC_CONTENT) {
    read_level_data_[index].assign({ScaleColor(GcContent(read), MaxGcContent)});
  } else if (channel_enum == DeepVariantChannelEnum::CH_IS_HOMOPOLYMER) {
    std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(read);
    read_level_data_[index] =
        ScaleColorVector(is_homopolymer, MaxIsHomoPolymer);
  } else if (channel_enum == DeepVariantChannelEnum::CH_HOMOPOLYMER_WEIGHTED) {
    std::vector<std::uint8_t> homopolymer_weighted = HomoPolymerWeighted(read);
    read_level_data_[index] =
        ScaleColorVector(homopolymer_weighted, MaxHomoPolymerWeighted);
  } else if (channel_enum == DeepVariantChannelEnum::CH_BLANK) {
    read_level_data_[index] = Blank(read);
  } else if (channel_enum == DeepVariantChannelEnum::CH_INSERT_SIZE) {
    read_level_data_[index] = ReadInsertSize(read);
  }

  return true;
}

// Calculate values for channels that depend on information at the
// granularity of bases within the read and/or reference sequence.
bool Channels::CalculateBaseLevelData(
    const Read& read, std::function<bool(int, int, const CigarUnit::Operation&)>
                          action_per_cigar_unit) {
  // In the following, we iterate over alignment information for each
  // base of read creating an association between a reference index,
  // read index and cigar operation and storing in a vector for subsequent
  // base level channel value calculation(s).
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
  bool ok = true;
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

    if (!ok) {
      return false;
    }
  }

  return true;
}

std::uint8_t Channels::GetChannelData(const std::string& channel, int col) {
  DeepVariantChannelEnum channel_enum = ChannelStrToEnum(channel);
  int index = channel_enum_to_index_[channel_enum];
  return data_[index][col];
}

void Channels::CalculateRefRows(
    const std::vector<DeepVariantChannelEnum>& channel_enums,
    const std::string& ref_bases) {
  int maxEnumValue = getMaxEnumValue(channel_enums);
  channel_enum_to_index_ = std::vector<int>(maxEnumValue + 1);
  int currIndex = 0;
  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    channel_enum_to_index_[channel_enum] = currIndex;
    currIndex++;
  }
  ref_data_ = std::vector<std::vector<unsigned char>>(channel_enums.size());
  // Calculates reference row values for each channel
  // Create a fake read to represent reference bases.
  Read refRead;
  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    int index = channel_enum_to_index_[channel_enum];
    if (channel_enum == DeepVariantChannelEnum::CH_READ_BASE) {
      ref_data_[index] = BaseColorVector(ref_bases, options_);
    } else if (channel_enum == DeepVariantChannelEnum::CH_BASE_QUALITY) {
      int ref_qual = options_.reference_base_quality();
      ref_data_[index].assign(
          ref_bases.size(), ScaleColor(ref_qual, options_.base_quality_cap()));
    } else if (channel_enum == DeepVariantChannelEnum::CH_MAPPING_QUALITY) {
      int ref_qual = options_.reference_base_quality();
      ref_data_[index].assign(
          ref_bases.size(), ScaleColor(ref_qual, options_.base_quality_cap()));
    } else if (channel_enum == DeepVariantChannelEnum::CH_STRAND) {
      int strand = StrandColor(true, options_);
      ref_data_[index].assign(ref_bases.size(),
                              static_cast<std::uint8_t>(strand));
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT) {
      int alt = SupportsAltColor(0, options_);
      ref_data_[index].assign(ref_bases.size(), static_cast<std::uint8_t>(alt));
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF) {
      int ref = MatchesRefColor(true, options_);
      ref_data_[index].assign(ref_bases.size(), static_cast<std::uint8_t>(ref));
    } else if (channel_enum == DeepVariantChannelEnum::CH_HAPLOTYPE_TAG) {
      ref_data_[index].assign(ref_bases.size(), {ScaleColor(0, 2)});
    } else if (channel_enum == DeepVariantChannelEnum::CH_ALLELE_FREQUENCY) {
      int allele_frequency_color = AlleleFrequencyColor(0, options_);
      ref_data_[index].assign(
          ref_bases.size(),
          {static_cast<std::uint8_t>(allele_frequency_color)});
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT) {
      ref_data_[index].assign(ref_bases.size(),
                              static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
    } else if (channel_enum == DeepVariantChannelEnum::CH_AVG_BASE_QUALITY) {
      ref_data_[index].assign(
          {static_cast<std::uint8_t>(kMaxPixelValueAsFloat)});
    } else if (channel_enum == DeepVariantChannelEnum::CH_IDENTITY) {
      ref_data_[index].assign(ref_bases.size(),
                              static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY) {
      ref_data_[index].assign(ref_bases.size(),
                              static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
    } else if (channel_enum == DeepVariantChannelEnum::CH_INSERT_SIZE) {
      ref_data_[index].assign(ref_bases.size(),
                              static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
    } else if (channel_enum == DeepVariantChannelEnum::CH_GC_CONTENT) {
      refRead.set_aligned_sequence(ref_bases);
      ref_data_[index].assign(ref_bases.size(),
                              ScaleColor(GcContent(refRead), MaxGcContent));
    } else if (channel_enum == DeepVariantChannelEnum::CH_IS_HOMOPOLYMER) {
      refRead.set_aligned_sequence(ref_bases);
      std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(refRead);
      ref_data_[index] = ScaleColorVector(is_homopolymer, MaxIsHomoPolymer);
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_HOMOPOLYMER_WEIGHTED) {
      refRead.set_aligned_sequence(ref_bases);
      std::vector<std::uint8_t> homopolymer_weighted =
          HomoPolymerWeighted(refRead);
      ref_data_[index] =
          ScaleColorVector(homopolymer_weighted, MaxHomoPolymerWeighted);
    } else {
      ref_data_[index].assign(ref_bases.size(), 0);
    }
  }
}

std::uint8_t Channels::GetRefRows(DeepVariantChannelEnum channel_enum,
                                  int col) {
  // Returns first value if size 1 else return specific column.
  // Note that ref_data is indexed by col and not pos.
  int index = channel_enum_to_index_[channel_enum];
  if (ref_data_[index].size() == 1) {
    return ref_data_[index][0];
  } else {
    return ref_data_[index][col];
  }
}

int Channels::getMaxEnumValue(
    const std::vector<DeepVariantChannelEnum>& channel_enums) {
  // Get the max channel enum value we need to support, used to size vectors
  // returns -1 if channel_enums is empty
  int maxEnumValue = -1;
  for (int i = 0; i < channel_enums.size(); ++i) {
    if (channel_enums[i] > maxEnumValue) {
      maxEnumValue = channel_enums[i];
    }
  }
  return maxEnumValue;
}
// Scales an input value to pixel range 0-254.
std::uint8_t Channels::ScaleColor(int value, float max_val) {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}

// Scales an input vector to pixel range 0-254
std::vector<std::uint8_t> Channels::ScaleColorVector(
    std::vector<std::uint8_t>& channel_values, float max_val) {
  for (int i = 0; i < channel_values.size(); i++) {
    int value = channel_values[i];
    if (static_cast<float>(value) > max_val) {
      value = max_val;
    }
    channel_values[i] = static_cast<int>(kMaxPixelValueAsFloat *
                                         (static_cast<float>(value) / max_val));
  }
  return channel_values;
}
//---------------//
// Base Channels //
//---------------//
int Channels::BaseColor(char base, const PileupImageOptions& options) {
  switch (base) {
    case 'A':
      return (options.base_color_offset_a_and_g() +
              options.base_color_stride() * 3);
    case 'G':
      return (options.base_color_offset_a_and_g() +
              options.base_color_stride() * 2);
    case 'T':
      return (options.base_color_offset_t_and_c() +
              options.base_color_stride() * 1);
    case 'C':
      return (options.base_color_offset_t_and_c() +
              options.base_color_stride() * 0);
    default:
      return 0;
  }
}

std::vector<std::uint8_t> Channels::BaseColorVector(
    const std::string& bases, const PileupImageOptions& options) {
  std::vector<std::uint8_t> base_colors;
  base_colors.reserve(bases.size());
  for (const char base : bases) {
    int color = BaseColor(base, options);
    base_colors.push_back(color);
  }
  return base_colors;
}

int Channels::StrandColor(bool on_positive_strand,
                          const PileupImageOptions& options) {
  return (on_positive_strand ? options.positive_strand_color()
                             : options.negative_strand_color());
}

int Channels::SupportsAltColor(int read_supports_alt,
                               const PileupImageOptions& options) {
  float alpha;
  if (read_supports_alt == 0) {
    alpha = options.allele_unsupporting_read_alpha();
  } else if (read_supports_alt == 1) {
    alpha = options.allele_supporting_read_alpha();
  } else {
    CHECK_EQ(read_supports_alt, 2) << "read_supports_alt can only be 0/1/2.";
    alpha = options.other_allele_supporting_read_alpha();
  }
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

// Does this read support ref, one of the alternative alleles, or an allele we
// aren't considering?
int Channels::ReadSupportsAlt(const DeepVariantCall& dv_call, const Read& read,
                              const std::vector<std::string>& alt_alleles) {
  std::string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
  for (const std::string& alt_allele : dv_call.variant().alternate_bases()) {
    const auto& allele_support = dv_call.allele_support();
    const bool alt_allele_present_in_call =
        allele_support.find(alt_allele) != allele_support.cend();

    if (alt_allele_present_in_call) {
      const auto& supp_read_names = allele_support.at(alt_allele).read_names();
      for (const std::string& read_name : supp_read_names) {
        const bool alt_in_alt_alleles =
            std::find(alt_alleles.begin(), alt_alleles.end(), alt_allele) !=
            alt_alleles.end();
        // Read can support an alt we are currently considering (1), a different
        // alt not present in alt_alleles (2), or ref (0).
        if (read_name == key && alt_in_alt_alleles) {
          return 1;
        } else if (read_name == key && !alt_in_alt_alleles) {
          return 2;
        }
      }
    }
  }
  return 0;
}

// Returns a value based on whether the current read base matched the
// reference base it was compared to.
int Channels::MatchesRefColor(bool base_matches_ref,
                              const PileupImageOptions& options) {
  float alpha = (base_matches_ref ? options.reference_matching_read_alpha()
                                  : options.reference_mismatching_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

//-----------------------//
// Experimental Channels //
//-----------------------//
constexpr bool Channels::isBaseLevelChannel(
    DeepVariantChannelEnum channelEnum) {
  switch (channelEnum) {
    case DeepVariantChannelEnum::CH_READ_BASE:
    case DeepVariantChannelEnum::CH_BASE_QUALITY:
    case DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF:
      return true;
    default:
      return false;
  }
}

DeepVariantChannelEnum Channels::ChannelStrToEnum(const std::string& channel) {
  if (channel == ch_read_base) return DeepVariantChannelEnum::CH_READ_BASE;
  if (channel == ch_base_quality)
    return DeepVariantChannelEnum::CH_BASE_QUALITY;
  if (channel == ch_mapping_quality)
    return DeepVariantChannelEnum::CH_MAPPING_QUALITY;
  if (channel == ch_strand) return DeepVariantChannelEnum::CH_STRAND;
  if (channel == ch_read_supports_variant)
    return DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT;
  if (channel == ch_base_differs_from_ref)
    return DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF;
  if (channel == ch_read_mapping_percent)
    return DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT;
  if (channel == ch_haplotype_tag)
    return DeepVariantChannelEnum::CH_HAPLOTYPE_TAG;
  if (channel == ch_allele_frequency)
    return DeepVariantChannelEnum::CH_ALLELE_FREQUENCY;
  if (channel == ch_diff_channels_alternate_allele_1)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_diff_channels_alternate_allele_2)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_avg_base_quality)
    return DeepVariantChannelEnum::CH_AVG_BASE_QUALITY;
  if (channel == ch_identity) return DeepVariantChannelEnum::CH_IDENTITY;
  if (channel == ch_gap_compressed_identity)
    return DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY;
  if (channel == ch_gc_content) return DeepVariantChannelEnum::CH_GC_CONTENT;
  if (channel == ch_is_homopolymer)
    return DeepVariantChannelEnum::CH_IS_HOMOPOLYMER;
  if (channel == ch_homopolymer_weighted)
    return DeepVariantChannelEnum::CH_HOMOPOLYMER_WEIGHTED;
  if (channel == ch_blank) return DeepVariantChannelEnum::CH_BLANK;
  if (channel == ch_insert_size) return DeepVariantChannelEnum::CH_INSERT_SIZE;
  if (channel == ch_base_channels_alternate_allele_1)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_base_channels_alternate_allele_2)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  CHECK(false) << "Channel '" << channel << "' should have a corresponding "
               << "enum in DeepVariantChannelEnum.";
}

std::vector<std::uint8_t> Channels::ReadInsertSize(const Read& read) {
  // Generates a vector reflecting the fragment length of the read
  std::vector<std::uint8_t> reads_with_insert_size(
      read.aligned_sequence().size(), normalizeFragmentLength(read));
  return reads_with_insert_size;
}

// normalizes a Read's `fragment_length` to a pixel value
int Channels::normalizeFragmentLength(const Read& read) {
  int fragment_length = std::abs(read.fragment_length());
  if (static_cast<float>(fragment_length) > MaxFragmentLength) {
    fragment_length = static_cast<int>(MaxFragmentLength);
  }
  return static_cast<int>(
      kMaxPixelValueAsFloat *
      (static_cast<float>(fragment_length) / MaxFragmentLength));
}

bool Channels::channel_exists(std::vector<std::string>& channels,
                              absl::string_view channel_name) {
  if (std::find(channels.begin(), channels.end(), channel_name) !=
      channels.end()) {
    return true;
  }
  return false;
}

int Channels::GetHPValueForHPChannel(const Read& read,
                                     int hp_tag_for_assembly_polishing) {
  // HP values are added to reads by DeepVariant (direct phasing).
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
  // If hp_tag_for_assembly_polishing is set to 2, this is a special case
  // assembly polishing:
  // If we're calling HP=2, when displayed with `--channel_list=haplotype`,
  // we want to swap the color of reads with HP=2 and HP=1.
  if (hp_tag_for_assembly_polishing == 2) {
    if (hp_value == 1) return 2;
    if (hp_value == 2) return 1;
  }

  // Otherwise, keep the default behavior.
  return hp_value;
}

// Get allele frequency color for a read.
// Convert a frequency value in float to color intensity (int) and normalize.
unsigned char Channels::AlleleFrequencyColor(
    float allele_frequency, const PileupImageOptions& options) {
  if (allele_frequency <= options.min_non_zero_allele_frequency()) {
    return 0;
  } else {
    float log10_af = log10(allele_frequency);
    float log10_min = log10(options.min_non_zero_allele_frequency());
    return ((log10_min - log10_af) / log10_min) *
           static_cast<int>(kMaxPixelValueAsFloat);
  }
}

// Get the allele frequency of the alt allele that is carried by a read.
float Channels::ReadAlleleFrequency(
    const DeepVariantCall& dv_call, const Read& read,
    const std::vector<std::string>& alt_alleles) {
  std::string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
  for (const std::string& alt_allele : dv_call.variant().alternate_bases()) {
    const auto& allele_support = dv_call.allele_support();
    auto it_read = allele_support.find(alt_allele);

    if (it_read != allele_support.end()) {
      const auto& supp_read_names = it_read->second.read_names();
      for (const std::string& read_name : supp_read_names) {
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

std::vector<std::uint8_t> Channels::Blank(const Read& read) {
  // Used to return a blank channel.
  std::vector<std::uint8_t> blank(read.aligned_sequence().size(), 0);
  return blank;
}

std::vector<std::uint8_t> Channels::HomoPolymerWeighted(const Read& read) {
  // Generates a vector reflecting the number of repeats observed.
  // ATCGGGAA
  // 11133322
  std::vector<std::uint8_t> homopolymer(read.aligned_sequence().size());
  auto seq = read.aligned_sequence();
  homopolymer[0] = 1;
  int current_weight = 1;
  for (int i = 1; i <= seq.size(); i++) {
    if (seq[i] == seq[i - 1]) {
      current_weight += 1;
    } else {
      for (int cw = current_weight; cw >= 1; cw--) {
        homopolymer[i - cw] = current_weight;
      }
      current_weight = 1;
    }
  }
  return homopolymer;
}

std::vector<std::uint8_t> Channels::IsHomoPolymer(const Read& read) {
  // Generates a vector indicating homopolymers of 3 or more.
  // ATCGGGAG
  // 00011100
  std::vector<std::uint8_t> homopolymer(read.aligned_sequence().size());
  auto seq = read.aligned_sequence();
  for (int i = 2; i < seq.size(); i++) {
    if (seq[i] == seq[i - 1] && seq[i - 1] == seq[i - 2]) {
      homopolymer[i] = 1;
      homopolymer[i - 1] = 1;
      homopolymer[i - 2] = 1;
    }
  }
  return homopolymer;
}

int Channels::GcContent(const Read& read) {
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

// Gap Compressed Identity: Ins/Del treated as individual events.
int Channels::GapCompressedIdentity(const Read& read) {
  // Calculates percentage of the read mapped to the reference
  int match_len = 0;
  int gap_compressed_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        gap_compressed_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        gap_compressed_len += op_len;
        break;
      case CigarUnit::INSERT:
        // Add a single event for insertion.
        gap_compressed_len += 1;
        break;
      case CigarUnit::DELETE:
        // Add a single event for a deletion.
        gap_compressed_len += 1;
        break;
      default:
        break;
    }
  }
  float gap_compressed_identity = static_cast<float>(match_len) /
                                  static_cast<float>(gap_compressed_len) * 100;
  return static_cast<int>(gap_compressed_identity);
}

// Identity: Similar to mapping percent but with a slightly different def.
int Channels::Identity(const Read& read) {
  int match_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        break;
      case CigarUnit::INSERT:
        break;
      case CigarUnit::DELETE:
        break;
      default:
        break;
    }
  }
  float mapping_percent = (static_cast<float>(match_len) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100;
  return static_cast<int>(mapping_percent);
}

// Average Base Quality: Averages base quality over length of read.
int Channels::AvgBaseQuality(const Read& read) {
  int base_qual_sum = 0;
  for (const auto& base_qual : read.aligned_quality()) {
    base_qual_sum += base_qual;
    // Base qualities range between 0 and 93
    if (base_qual < 0 || base_qual > 93) {
      LOG(FATAL) << "Encountered base quality outside of bounds (0,93):"
                 << base_qual << ", read=" << read.fragment_name();
    }
  }
  float avg_base_qual = (static_cast<float>(base_qual_sum) /
                         static_cast<float>(read.aligned_quality().size()));
  return static_cast<int>(avg_base_qual);
}

// Read Mapping Percent: Calculates percentage of bases mapped to reference.
int Channels::ReadMappingPercent(const Read& read) {
  int match_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        break;
      default:
        break;
    }
  }
  float mapping_percent = (static_cast<float>(match_len) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100;
  return static_cast<int>(mapping_percent);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
