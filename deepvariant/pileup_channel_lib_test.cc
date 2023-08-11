/*
 * Copyright 2021 Google LLC.
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

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

using nucleus::genomics::v1::Read;

namespace learning {
namespace genomics {
namespace deepvariant {

TEST(ScaleColor, BasicCase) {
  EXPECT_EQ(ScaleColor(50, 100), 127);
  EXPECT_EQ(ScaleColor(127, kMaxPixelValueAsFloat), 127);
  // Beyond max scales to max.
  EXPECT_EQ(ScaleColor(500, kMaxPixelValueAsFloat), kMaxPixelValueAsFloat);
}

TEST(ScaleColorVector, BasicCase) {
  std::vector<std::uint8_t> test_vector{5, 10, 25};
  std::vector<std::uint8_t> expect_vector{25, 50, 127};
  EXPECT_EQ(ScaleColorVector(test_vector, 50), expect_vector);
}

TEST(ScaleColorVectorLarge, OverMaxCase) {
  std::vector<std::uint8_t> test_vector;
  test_vector.resize(500);
  for (int i = 0; i < test_vector.size(); i++) {
    test_vector[i] = i + 1;
  }
  test_vector = ScaleColorVector(test_vector, kMaxPixelValueAsFloat);
  std::uint8_t j = 0;
  for (auto& i : test_vector) {
    j++;
    if (i < kMaxPixelValueAsFloat) {
      EXPECT_EQ(i, j);
    } else {
      EXPECT_EQ(i, static_cast<int>(kMaxPixelValueAsFloat));
    }
  }
}

TEST(BaseColor, A) {
  PileupImageOptions options{};
  options.set_base_color_offset_a_and_g(1);
  options.set_base_color_stride(1);
  std::uint8_t color = BaseColor('A', options);
  EXPECT_EQ(color, 4);
}

TEST(BaseColor, T) {
  PileupImageOptions options{};
  options.set_base_color_offset_t_and_c(1);
  options.set_base_color_stride(1);
  std::uint8_t color = BaseColor('T', options);
  EXPECT_EQ(color, 2);
}

TEST(BaseColor, G) {
  PileupImageOptions options{};
  options.set_base_color_offset_a_and_g(1);
  options.set_base_color_stride(1);
  std::uint8_t color = BaseColor('G', options);
  EXPECT_EQ(color, 3);
}

TEST(BaseColor, C) {
  PileupImageOptions options{};
  options.set_base_color_offset_t_and_c(1);
  options.set_base_color_stride(1);
  std::uint8_t color = BaseColor('C', options);
  EXPECT_EQ(color, 1);
}

TEST(BaseColorVector, ATGC) {
  PileupImageOptions options{};
  options.set_base_color_offset_a_and_g(1);
  options.set_base_color_offset_t_and_c(1);
  options.set_base_color_stride(1);

  std::vector<std::uint8_t> expect_vector{4, 2, 3, 1};
  std::vector<std::uint8_t> colors_vector = BaseColorVector("ATGC", options);
  EXPECT_EQ(colors_vector, expect_vector);
}

TEST(StrandColor, PositiveStrand) {
  PileupImageOptions options{};
  options.set_positive_strand_color(10);
  options.set_negative_strand_color(20);
  std::uint8_t sc = StrandColor(true, options);
  EXPECT_EQ(sc, 10);
}

TEST(StrandColor, NegativeStrand) {
  PileupImageOptions options{};
  options.set_positive_strand_color(10);
  options.set_negative_strand_color(20);
  std::uint8_t sc = StrandColor(false, options);
  EXPECT_EQ(sc, 20);
}

TEST(ReadSupportsAlt, AlleleUnsupporting) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTT", {"8M"});
  DeepVariantCall dv_call = DeepVariantCall::default_instance();
  std::vector<std::string> alt_alleles = {};

  std::uint8_t rsa = ReadSupportsAlt(dv_call, read, alt_alleles);
  EXPECT_EQ(rsa, 0);
}

TEST(ReadSupportsAlt, AlleleSupporting) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTT", {"8M"}, "FRAG1");
  read.set_read_number(1);

  DeepVariantCall_SupportingReads dv_supporting_reads =
      DeepVariantCall_SupportingReads::default_instance();
  dv_supporting_reads.add_read_names("FRAG1/1");

  DeepVariantCall dv_call = DeepVariantCall::default_instance();

  dv_call.mutable_variant()->mutable_alternate_bases()->Add("GGGCGCATT");

  dv_call.mutable_allele_support()->insert(
      google::protobuf::MapPair<std::string, DeepVariantCall_SupportingReads>(
          "GGGCGCATT", dv_supporting_reads));

  std::vector<std::string> alt_alleles = {"GGGCGCATT"};

  std::uint8_t rsa = ReadSupportsAlt(dv_call, read, alt_alleles);
  EXPECT_EQ(rsa, 1);
}

TEST(ReadSupportsAlt, OtherAlleleSupporting) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTT", {"8M"}, "FRAG2");
  read.set_read_number(2);

  DeepVariantCall_SupportingReads dv_supporting_reads =
      DeepVariantCall_SupportingReads::default_instance();
  dv_supporting_reads.add_read_names("FRAG2/2");

  DeepVariantCall dv_call = DeepVariantCall::default_instance();

  dv_call.mutable_variant()->mutable_alternate_bases()->Add("GGGCGCATT");

  dv_call.mutable_allele_support()->insert(
      google::protobuf::MapPair<std::string, DeepVariantCall_SupportingReads>(
          "GGGCGCATT", dv_supporting_reads));

  std::vector<std::string> alt_alleles = {};

  std::uint8_t rsa = ReadSupportsAlt(dv_call, read, alt_alleles);
  EXPECT_EQ(rsa, 2);
}

TEST(MatchesRefColor, BaseMatch) {
  PileupImageOptions options{};
  options.set_reference_matching_read_alpha(1);
  options.set_reference_mismatching_read_alpha(0);
  int match = MatchesRefColor(true, options);
  EXPECT_EQ(match, kMaxPixelValueAsFloat);
}

TEST(MatchesRefColor, BaseMistmatch) {
  PileupImageOptions options{};
  options.set_reference_matching_read_alpha(0);
  options.set_reference_mismatching_read_alpha(1);
  int match = MatchesRefColor(false, options);
  EXPECT_EQ(match, kMaxPixelValueAsFloat);
}

TEST(SupportsAltColor, AlleleSupporting) {
  PileupImageOptions options{};
  options.set_allele_unsupporting_read_alpha(1.0);
  options.set_allele_supporting_read_alpha(0.0);
  options.set_other_allele_supporting_read_alpha(0.0);
  std::uint8_t sac = SupportsAltColor(0, options);
  EXPECT_EQ(sac, kMaxPixelValueAsFloat);
}

TEST(SupportsAltColor, AlleleUnsupporting) {
  PileupImageOptions options{};
  options.set_allele_unsupporting_read_alpha(0.0);
  options.set_allele_supporting_read_alpha(1.0);
  options.set_other_allele_supporting_read_alpha(0.0);
  std::uint8_t sac = SupportsAltColor(1, options);
  EXPECT_EQ(sac, kMaxPixelValueAsFloat);
}

TEST(SupportsAltColor, OtherAlleleSupporting) {
  PileupImageOptions options{};
  options.set_allele_unsupporting_read_alpha(0.0);
  options.set_allele_supporting_read_alpha(0.0);
  options.set_other_allele_supporting_read_alpha(1.0);
  std::uint8_t sac = SupportsAltColor(2, options);
  EXPECT_EQ(sac, kMaxPixelValueAsFloat);
}

TEST(ReadMappingPercentTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5M", "5D"});
  std::uint8_t rmp = ReadMappingPercent(read);
  EXPECT_EQ(rmp, 50);
}

TEST(AvgBaseQualityTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"10M"});
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, i + 1);
  }
  std::uint8_t rmp = AvgBaseQuality(read);
  EXPECT_EQ(rmp, 5);
}

TEST(AvgBaseQualityTest, BaseQualityTooHigh) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"10M"});
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, 100);
  }
  EXPECT_DEATH(AvgBaseQuality(read),
               "Encountered base quality outside of bounds");
}

TEST(IdentityTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5M", "1I", "4M"});
  std::uint8_t id = Identity(read);
  EXPECT_EQ(id, 90);
}

TEST(IdentityTest, PacBioStyleCigar) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5=", "1X", "4="});
  std::uint8_t id = Identity(read);
  EXPECT_EQ(id, 90);
}

TEST(GapCompressedIdentityTest, InsertionCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3M", "4I", "3M"});
  std::uint8_t id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 85);
}

TEST(GapCompressedIdentityTest, DeletionCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3M", "4D", "3M"});
  std::uint8_t id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 85);
}

TEST(GapCompressedIdentityTest, PacBioStyleCigar) {
  Read read =
      nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3=", "2X", "2I", "3="});
  std::uint8_t id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 66);
}

TEST(GcContestTest, AllGc) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGGGCCCCC", {"10M"});
  std::uint8_t gc = GcContent(read);
  EXPECT_EQ(gc, 100);
}

TEST(GcContestTest, HalfGc) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGGGTTTTT", {"10M"});
  std::uint8_t gc = GcContent(read);
  EXPECT_EQ(gc, 50);
}

TEST(IsHomoPolymerTest, IsHomopolymerBeginning) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGATAATA", {"9M"});
  std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(read);
  std::vector<std::uint8_t> expected{1, 1, 1, 0, 0, 0, 0, 0, 0};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerMiddle) {
  Read read = nucleus::MakeRead("chr1", 1, "ATTGGGTTA", {"9M"});
  std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(read);
  std::vector<std::uint8_t> expected{0, 0, 0, 1, 1, 1, 0, 0, 0};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerEnd) {
  Read read = nucleus::MakeRead("chr1", 1, "ATAATAGGG", {"9M"});
  std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(read);
  std::vector<std::uint8_t> expected{0, 0, 0, 0, 0, 0, 1, 1, 1};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerAll) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAAAAAA", {"9M"});
  std::vector<std::uint8_t> is_homopolymer = IsHomoPolymer(read);
  std::vector<std::uint8_t> expected{1, 1, 1, 1, 1, 1, 1, 1, 1};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(HomoPolymerWeightedTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  std::vector<std::uint8_t> w_homopolymer = HomoPolymerWeighted(read);
  std::vector<std::uint8_t> expected{1, 1, 2, 2, 3, 3, 3, 4,
                                     4, 4, 4, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_homopolymer, expected);
}

TEST(HomoPolymerWeightedTest, WeightedHomoPolymerMax) {
  std::string bases;
  bases.insert(0, 20, 'A');
  bases.insert(0, 10, 'G');
  Read read = nucleus::MakeRead("chr1", 1, bases, {"60M"});
  std::vector<std::uint8_t> w_homopolymer = HomoPolymerWeighted(read);
  std::vector<std::uint8_t> expected;
  expected.insert(expected.begin(), 20, 20);
  expected.insert(expected.begin(), 10, 10);
  EXPECT_EQ(w_homopolymer, expected);
}

TEST(ReadInsertSizeTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(22);
  std::vector<std::uint8_t> w_insert_size = ReadInsertSize(read);
  std::vector<std::uint8_t> expected{5, 5, 5, 5, 5, 5, 5, 5,
                                     5, 5, 5, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, NegativeValueCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(-22);
  std::vector<std::uint8_t> w_insert_size = ReadInsertSize(read);
  std::vector<std::uint8_t> expected{5, 5, 5, 5, 5, 5, 5, 5,
                                     5, 5, 5, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, ExceedsLimit) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(1001);
  std::vector<std::uint8_t> w_insert_size = ReadInsertSize(read);
  std::vector<std::uint8_t> expected{254, 254, 254, 254, 254, 254, 254, 254,
                                     254, 254, 254, 254, 254, 254, 254, 254};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, NoValue) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  std::vector<std::uint8_t> w_insert_size = ReadInsertSize(read);
  std::vector<std::uint8_t> expected{0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(GetChannelDataTest, ReadData) {
  PileupImageOptions options{};
  options.set_mapping_quality_cap(1);
  options.set_positive_strand_color(20);
  options.set_allele_unsupporting_read_alpha(1.0);
  options.set_base_color_offset_a_and_g(1);
  options.set_base_color_offset_t_and_c(1);
  options.set_base_color_stride(1);
  options.set_base_quality_cap(20);
  options.set_reference_matching_read_alpha(1);
  options.set_reference_mismatching_read_alpha(0);

  OptChannels channel_set{options};
  std::vector<std::string> channels{ch_read_base,
                                    ch_base_quality,
                                    ch_mapping_quality,
                                    ch_strand,
                                    ch_read_supports_variant,
                                    ch_base_differs_from_ref,
                                    ch_read_mapping_percent,
                                    ch_avg_base_quality,
                                    ch_identity,
                                    ch_gap_compressed_identity,
                                    ch_gc_content,
                                    ch_is_homopolymer,
                                    ch_homopolymer_weighted,
                                    ch_blank,
                                    ch_insert_size};
  Read ref_read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTTAT", {"11M"});
  Read read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTTAT", {"11M"});
  read.set_fragment_length(MaxFragmentLength);
  const int base_quality = 33;
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, base_quality);
  }

  DeepVariantCall dv_call = DeepVariantCall::default_instance();
  std::vector<std::string> alt_alleles = {};

  channel_set.CalculateChannels(channels, read, ref_read.aligned_sequence(),
                                dv_call, alt_alleles, 0);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_base, 11), 4);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_base, 9), 2);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_base, 1), 3);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_base, 4), 1);
  EXPECT_EQ(channel_set.GetChannelData(ch_base_quality, 1),
            kMaxPixelValueAsFloat);
  EXPECT_EQ(channel_set.GetChannelData(ch_mapping_quality, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_strand, 1), 20);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_supports_variant, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_base_differs_from_ref, 1),
            kMaxPixelValueAsFloat);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_mapping_percent, 3), 231);
  EXPECT_EQ(channel_set.GetChannelData(ch_avg_base_quality, 3), 90);
  EXPECT_EQ(channel_set.GetChannelData(ch_identity, 9), 231);
  EXPECT_EQ(channel_set.GetChannelData(ch_gap_compressed_identity, 9),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_gc_content, 3), 127);
  EXPECT_EQ(channel_set.GetChannelData(ch_is_homopolymer, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_is_homopolymer, 4), 0);
  EXPECT_EQ(channel_set.GetChannelData(ch_homopolymer_weighted, 1), 25);
  EXPECT_EQ(channel_set.GetChannelData(ch_homopolymer_weighted, 9), 33);
  EXPECT_EQ(channel_set.GetChannelData(ch_blank, 1), 0);
  EXPECT_EQ(channel_set.GetChannelData(ch_insert_size, 1), 254);
}

TEST(GetRefChannelDataTest, ReadData) {
  PileupImageOptions options{};
  options.set_reference_base_quality(20);
  options.set_base_quality_cap(20);
  options.set_allele_unsupporting_read_alpha(1.0);
  options.set_positive_strand_color(20);
  options.set_base_color_offset_a_and_g(1);
  options.set_base_color_offset_t_and_c(1);
  options.set_base_color_stride(1);
  options.set_reference_matching_read_alpha(1.0);

  OptChannels channel_set{options};
  std::vector<std::string> channels{ch_read_base,
                                    ch_base_quality,
                                    ch_mapping_quality,
                                    ch_strand,
                                    ch_read_supports_variant,
                                    ch_base_differs_from_ref,
                                    ch_read_mapping_percent,
                                    ch_avg_base_quality,
                                    ch_identity,
                                    ch_gap_compressed_identity,
                                    ch_gc_content,
                                    ch_is_homopolymer,
                                    ch_homopolymer_weighted,
                                    ch_blank};
  Read ref_read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTTAT", {"11M"});
  const int base_quality = 33;
  for (size_t i = 0; i < ref_read.aligned_sequence().size(); ++i) {
    ref_read.set_aligned_quality(i, base_quality);
  }
  channel_set.CalculateRefRows(channels, ref_read.aligned_sequence());
  EXPECT_EQ(channel_set.GetRefRows(ch_read_base, 10), 4);
  EXPECT_EQ(channel_set.GetRefRows(ch_read_base, 8), 2);
  EXPECT_EQ(channel_set.GetRefRows(ch_read_base, 0), 3);
  EXPECT_EQ(channel_set.GetRefRows(ch_read_base, 3), 1);
  EXPECT_EQ(channel_set.GetRefRows(ch_base_quality, 0), kMaxPixelValueAsFloat);
  EXPECT_EQ(channel_set.GetRefRows(ch_mapping_quality, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_strand, 1), 20);
  EXPECT_EQ(channel_set.GetRefRows(ch_read_supports_variant, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_base_differs_from_ref, 0),
            kMaxPixelValueAsFloat);
  EXPECT_EQ(channel_set.GetRefRows(ch_read_mapping_percent, 3),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_avg_base_quality, 3),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_identity, 9), 254);
  EXPECT_EQ(channel_set.GetRefRows(ch_gap_compressed_identity, 9),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_gc_content, 3), 127);
  EXPECT_EQ(channel_set.GetRefRows(ch_is_homopolymer, 1),
            static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_is_homopolymer, 4), 0);
  EXPECT_EQ(channel_set.GetRefRows(ch_homopolymer_weighted, 1), 25);
  EXPECT_EQ(channel_set.GetRefRows(ch_homopolymer_weighted, 9), 33);
  EXPECT_EQ(channel_set.GetRefRows(ch_blank, 0), 0);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
