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

#include "deepvariant/pileup_image_native.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/testing_utils.h"
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Read = nucleus::genomics::v1::Read;

struct GetSampleAltAlignedPileupTestData {
  std::string test_name;
  AltAlignedPileup global_alt_aligned_pileup;
  std::string sample_alt_aligned_pileup_name;
  AltAlignedPileup expected_sample_alt_aligned_pileup;
};

class GetSampleAltAlignedPileupTest
    : public testing::TestWithParam<GetSampleAltAlignedPileupTestData> {};

TEST_P(GetSampleAltAlignedPileupTest, GetSampleAltAlignedPileupTestCases) {
  const GetSampleAltAlignedPileupTestData& param = GetParam();
  CHECK_EQ(GetSampleAltAlignedPileup(param.global_alt_aligned_pileup,
                                     param.sample_alt_aligned_pileup_name),
           param.expected_sample_alt_aligned_pileup);
}

INSTANTIATE_TEST_SUITE_P(
    GetSampleAltAlignedPileupTests, GetSampleAltAlignedPileupTest,
    testing::ValuesIn(std::vector<GetSampleAltAlignedPileupTestData>({
        {
            .test_name =
                "chooses_global_alt_aligned_pileup_when_sample_is_empty",
            .global_alt_aligned_pileup = AltAlignedPileup::kRows,
            .sample_alt_aligned_pileup_name = "",
            .expected_sample_alt_aligned_pileup = AltAlignedPileup::kRows,
        },
        {
            .test_name = "chooses_none_when_sample_specifies_none",
            .global_alt_aligned_pileup = AltAlignedPileup::kRows,
            .sample_alt_aligned_pileup_name = "none",
            .expected_sample_alt_aligned_pileup = AltAlignedPileup::kNone,
        },
        {
            .test_name = "chooses_single_row_when_sample_specifies_single_row",
            .global_alt_aligned_pileup = AltAlignedPileup::kRows,
            .sample_alt_aligned_pileup_name = "single_row",
            .expected_sample_alt_aligned_pileup = AltAlignedPileup::kSingleRow,
        },
        {
            .test_name = "global_is_equal_to_sample",
            .global_alt_aligned_pileup = AltAlignedPileup::kDiffChannels,
            .sample_alt_aligned_pileup_name = "diff_channels",
            .expected_sample_alt_aligned_pileup =
                AltAlignedPileup::kDiffChannels,
        },
        {
            .test_name = "global_is_none_and_sample_is_rows",
            .global_alt_aligned_pileup = AltAlignedPileup::kNone,
            .sample_alt_aligned_pileup_name = "rows",
            .expected_sample_alt_aligned_pileup = AltAlignedPileup::kRows,
        },
    })));

struct CalculatePileupImageHeightTestData {
  std::string test_name;
  std::string pic_options_alt_aligned_pileup;
  int sample_1_pileup_height;
  std::string sample_1_alt_aligned_pileup;
  int sample_2_pileup_height;
  std::string sample_2_alt_aligned_pileup;
  int expected_pileup_image_height;
};

class CalculatePileupImageHeightTest
    : public testing::TestWithParam<CalculatePileupImageHeightTestData> {};

TEST_P(CalculatePileupImageHeightTest, CalculatePileupImageHeightTestCases) {
  const CalculatePileupImageHeightTestData& param = GetParam();
  MakeExamplesOptions options;
  options.mutable_pic_options()->set_alt_aligned_pileup(
      param.pic_options_alt_aligned_pileup);

  // Add first sample
  SampleOptions* sample_options = options.add_sample_options();
  sample_options->set_name("sample_1");
  sample_options->set_pileup_height(param.sample_1_pileup_height);
  sample_options->set_alt_aligned_pileup(param.sample_1_alt_aligned_pileup);

  // Add second sample
  SampleOptions* sample_options_2 = options.add_sample_options();
  sample_options_2->set_name("sample_2");
  sample_options_2->set_pileup_height(param.sample_2_pileup_height);
  sample_options_2->set_alt_aligned_pileup(param.sample_2_alt_aligned_pileup);

  CHECK_EQ(CalculatePileupImageHeight(options),
           param.expected_pileup_image_height);
}

INSTANTIATE_TEST_SUITE_P(
    CalculatePileupImageHeightTests, CalculatePileupImageHeightTest,
    testing::ValuesIn(std::vector<CalculatePileupImageHeightTestData>({
        {
            .test_name = "no_alt_aligned_pileup_set",
            .pic_options_alt_aligned_pileup = "none",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "none",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "none",
            .expected_pileup_image_height = 20,
        },
        {
            .test_name = "sample_1_single_row",
            .pic_options_alt_aligned_pileup = "none",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "none",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "single_row",
            .expected_pileup_image_height = 30,
        },
        {
            .test_name = "global_rows_no_sample_alt_aligned_pileup_set",
            .pic_options_alt_aligned_pileup = "rows",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "",
            .expected_pileup_image_height = 60,
        },
        {
            .test_name = "global_rows_one_sample_alt_aligned_pileup_sets_none",
            .pic_options_alt_aligned_pileup = "rows",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "none",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "",
            .expected_pileup_image_height = 40,
        },
        {
            .test_name = "global_none_one_sample_alt_aligned_pileup_sets_rows",
            .pic_options_alt_aligned_pileup = "none",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "none",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "rows",
            .expected_pileup_image_height = 40,
        },
        {
            .test_name =
                "global_single_row_one_sample_alt_aligned_pileup_sets_none",
            .pic_options_alt_aligned_pileup = "single_row",
            .sample_1_pileup_height = 10,
            .sample_1_alt_aligned_pileup = "none",
            .sample_2_pileup_height = 10,
            .sample_2_alt_aligned_pileup = "",
            .expected_pileup_image_height = 30,
        },
    })));

struct BuildPileupForOneSampleTestData {
  std::string test_name;
  DeepVariantCall dv_call;
  std::string ref_bases;
  std::vector<Read> reads;
  int image_start_pos;
  std::vector<std::string> alt_alleles;
  std::vector<ImageRow> expected_image_rows;
};

class BuildPileupForOneSampleTest
    : public testing::TestWithParam<BuildPileupForOneSampleTestData> {};

TEST_P(BuildPileupForOneSampleTest, BuildPileupForOneSampleTestCases) {
  const BuildPileupForOneSampleTestData& param = GetParam();
  PileupImageOptions pileup_image_options = MakeDefaultPileupImageOptions(
      /*width=*/ 11, /*height=*/ 4, /*ref_band_height=*/ 1);
  auto channels = pileup_image_options.mutable_channels();
  channels->Add("read_base");
  channels->Add("base_quality");
  channels->Add("mapping_quality");
  pileup_image_options.set_num_channels(channels->size());

  std::vector<const Read*> read_ptrs;
  for (const Read& read : param.reads) {
    read_ptrs.push_back(&read);
  }

  std::vector<std::unique_ptr<ImageRow>> image_rows =
      PileupImageEncoderNative(pileup_image_options).BuildPileupForOneSample(
          param.dv_call,
          param.ref_bases,
          read_ptrs,
          param.image_start_pos,
          param.alt_alleles,
          SampleOptions(),
          0.0, nullptr, {});
  // Convert to a vector of ImageRow to make it easier to compare with the
  // expected image rows.
  std::vector<ImageRow> image_rows_vec;
  for (const std::unique_ptr<ImageRow>& image_row : image_rows) {
    image_rows_vec.push_back(*image_row);
  }

  // Output image rows in human readable format. This is useful for debugging
  // failed tests.
  for (const ImageRow& image_row : image_rows_vec) {
    for (int channel = 0; channel < image_row.num_channels; ++channel) {
      std::string channel_row;
      for (int pos = 0; pos < image_row.width; ++pos) {
        channel_row += ", " +
            std::to_string(image_row.channel_data[channel][pos]);
      }
      LOG(INFO) << channel_row;
    }
  }

  EXPECT_THAT(image_rows_vec,
              testing::UnorderedElementsAreArray(param.expected_image_rows));
};

// For all test cases the pileup image options are the same. Width is 11, height
// is 4, ref_band_height is 1.
// Tests compare the image rows in the same order as in the expected_image_rows.
// expected image rows are created using MakeImageRow, the first row is the
// reference band, the rest are the main part of the image.
INSTANTIATE_TEST_SUITE_P(
    BuildPileupForOneSampleTests, BuildPileupForOneSampleTest,
    testing::ValuesIn(std::vector<BuildPileupForOneSampleTestData>({
      {
        .test_name = "simple_case",
        .dv_call = MakeDeepVariantCall(MakeVariant("A", {"G"}, 5)),
        .ref_bases = "ACGTACTCCCA",
        .reads = {
          MakeRead("chr1", 0, "ACGTGCTCCCA", {"11M"}, "read_2"),
          MakeRead("chr1", 0, "ACGTGCTCCCA", {"11M"}, "read_3"),
        },
        .image_start_pos = 0,
        .alt_alleles = {"G"},
        .expected_image_rows = {
                *MakeImageRow({
                    {250, 30, 180, 100, 250, 30, 100, 30, 30, 30, 250},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 180, 30, 100, 30, 30, 30, 250},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 180, 30, 100, 30, 30, 30, 250},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
                              11, 3).get(),
        }
      }
      , {
        .test_name = "no_reads",
        .dv_call = MakeDeepVariantCall(MakeVariant("A", {"G"}, 5)),
        .ref_bases = "ACGTACTCCCA",
        .reads = {},
        .image_start_pos = 0,
        .alt_alleles = {"G"},
        .expected_image_rows = {
                *MakeImageRow({
                    {250, 30, 180, 100, 250, 30, 100, 30, 30, 30, 250},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
                              11, 3).get(),
        }
      }
      , {
        .test_name = "numer_of_reads_greater_than_max_reads",
        .dv_call = MakeDeepVariantCall(MakeVariant("A", {"AGG"}, 5)),
        .ref_bases = "ACGTACTCCCA",
        .reads = {
          MakeRead("chr1", 0, "ACGTAGGCTCCCA", {"5M", "2I", "5M"}, "read_2"),
          MakeRead("chr1", 0, "ACGTAGGCTCCCA", {"5M", "2I", "5M"}, "read_3"),
          MakeRead("chr1", 0, "ACGTAGGGCTCCCA", {"5M", "3I", "5M"}, "read_4"),
          MakeRead("chr1", 0, "ACGTAGGGCTCCCA", {"5M", "3I", "5M"}, "read_5"),
        },
        .image_start_pos = 0,
        .alt_alleles = {"G"},
        .expected_image_rows = {
                *MakeImageRow({
                    {250, 30, 180, 100, 250, 30, 100, 30, 30, 30, 250},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
        }
      }
      , {
        .test_name = "image_creation_with_haplotype_sorting",
        .dv_call = MakeDeepVariantCall(MakeVariant("A", {"AGG", "AGGG"}, 5)),
        .ref_bases = "ACGTACTCCCA",
        .reads = {
          MakeRead("chr1", 0, "ACGTAGGCTCCCA", {"5M", "2I", "5M"}, "read_2", 2),
          MakeRead("chr1", 0, "TCGTAGGCTCCCA", {"5M", "2I", "5M"}, "read_3", 0),
          MakeRead("chr1", 0, "CCGTAGGGCTCCCA", {"5M", "3I", "5M"}, "read_4",
                   1),
        },
        .image_start_pos = 0,
        .alt_alleles = {"G"},
        .expected_image_rows = {
                *MakeImageRow({
                    {250, 30, 180, 100, 250, 30, 100, 30, 30, 30, 250},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254}},
                              11, 3).get(),
                *MakeImageRow({
                    {250, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {100, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
                *MakeImageRow({
                    {30, 30, 180, 100, 0, 30, 100, 30, 30, 30, 0},
                    {190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 0},
                    {254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 0}},
                              11, 3).get(),
        }
      }
})));

std::vector<unsigned char> ConcatExpectedPileupImages(
    std::vector<std::string> expected_pileups,
    std::unordered_map<std::string, std::vector<unsigned char>>
        expected_pileups_by_name) {
  std::vector<unsigned char> concatenated_vector;
  for (std::string name : expected_pileups) {
    std::vector<unsigned char> row = expected_pileups_by_name[name];
    std::copy(row.begin(), row.end(), std::back_inserter(concatenated_vector));
  }
  return concatenated_vector;
}

struct FillPileupArrayTestData {
  std::string test_name;
  int expected_num_rows;
  int expected_num_channels;
  AltAlignedPileup alt_aligned_pileup;
  std::vector<int> alt_image_row_indices;
  std::vector<std::string> expected_pileup_images;
};

class FillPileupArrayTest
    : public testing::TestWithParam<FillPileupArrayTestData> {};

TEST_P(FillPileupArrayTest, FillPileupArrayTestCases) {
  int num_channels = 7;
  int width = 5;
  std::vector<std::unique_ptr<ImageRow>> input;
  input.emplace_back(MakeImageRow({{11, 11, 11, 11, 11},
                                   {21, 21, 21, 21, 21},
                                   {31, 31, 31, 31, 31},
                                   {41, 41, 41, 41, 41},
                                   {51, 51, 51, 51, 51},
                                   {61, 61, 61, 61, 61},
                                   {71, 71, 71, 71, 71}},
                                  width, num_channels));
  input.emplace_back(MakeImageRow({{12, 12, 12, 12, 12},
                                   {22, 22, 22, 22, 22},
                                   {32, 32, 32, 32, 32},
                                   {42, 42, 42, 42, 42},
                                   {52, 52, 52, 52, 52},
                                   {62, 62, 62, 62, 62},
                                   {72, 72, 72, 72, 72}},
                                  width, num_channels));
  input.emplace_back(MakeImageRow({{13, 13, 13, 13, 13},
                                   {23, 23, 23, 23, 23},
                                   {33, 33, 33, 33, 33},
                                   {43, 43, 43, 43, 43},
                                   {53, 53, 53, 53, 53},
                                   {63, 63, 63, 63, 63},
                                   {73, 73, 73, 73, 73}},
                                  width, num_channels));
  // In order to distinguish channel 8 and 9 "-5" is added to channel 8 and
  // "+5" is added to channel 9.
  std::vector<std::vector<std::unique_ptr<ImageRow>>> alt_rows(2);
  alt_rows[0].emplace_back(
      MakeImageRow({{11 - 5, 11 - 5, 11 - 5, 11 - 5, 11 - 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {61 - 5, 61 - 5, 61 - 5, 61 - 5, 61 - 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));
  alt_rows[0].emplace_back(
      MakeImageRow({{12 - 5, 12 - 5, 12 - 5, 12 - 5, 12 - 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {62 - 5, 62 - 5, 62 - 5, 62 - 5, 62 - 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));
  alt_rows[0].emplace_back(
      MakeImageRow({{13 - 5, 13 - 5, 13 - 5, 13 - 5, 13 - 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {63 - 5, 63 - 5, 63 - 5, 63 - 5, 63 - 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));
  std::vector<std::unique_ptr<ImageRow>> alt_rows_2;
  alt_rows[1].emplace_back(
      MakeImageRow({{11 + 5, 11 + 5, 11 + 5, 11 + 5, 11 + 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {61 + 5, 61 + 5, 61 + 5, 61 + 5, 61 + 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));
  alt_rows[1].emplace_back(
      MakeImageRow({{12 + 5, 12 + 5, 12 + 5, 12 + 5, 12 + 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {62 + 5, 62 + 5, 62 + 5, 62 + 5, 62 + 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));
  alt_rows[1].emplace_back(
      MakeImageRow({{13 + 5, 13 + 5, 13 + 5, 13 + 5, 13 + 5},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1},
                    {63 + 5, 63 + 5, 63 + 5, 63 + 5, 63 + 5},
                    {1, 1, 1, 1, 1}},
                   width, num_channels));

  std::unordered_map<std::string, std::vector<unsigned char>>
      expected_pileup_image_by_name = {
          {"ref_row",
           {
               11, 21, 31, 41, 51, 61, 71, 11, 21, 31, 41, 51,
               61, 71, 11, 21, 31, 41, 51, 61, 71, 11, 21, 31,
               41, 51, 61, 71, 11, 21, 31, 41, 51, 61, 71,

               12, 22, 32, 42, 52, 62, 72, 12, 22, 32, 42, 52,
               62, 72, 12, 22, 32, 42, 52, 62, 72, 12, 22, 32,
               42, 52, 62, 72, 12, 22, 32, 42, 52, 62, 72,

               13, 23, 33, 43, 53, 63, 73, 13, 23, 33, 43, 53,
               63, 73, 13, 23, 33, 43, 53, 63, 73, 13, 23, 33,
               43, 53, 63, 73, 13, 23, 33, 43, 53, 63, 73,
           }},
          {"alt_row_1",
           {
               11 - 5, 1, 1,      1, 1,      61 - 5, 1, 11 - 5, 1, 1,      1, 1,
               61 - 5, 1, 11 - 5, 1, 1,      1,      1, 61 - 5, 1, 11 - 5, 1, 1,
               1,      1, 61 - 5, 1, 11 - 5, 1,      1, 1,      1, 61 - 5, 1,

               12 - 5, 1, 1,      1, 1,      62 - 5, 1, 12 - 5, 1, 1,      1, 1,
               62 - 5, 1, 12 - 5, 1, 1,      1,      1, 62 - 5, 1, 12 - 5, 1, 1,
               1,      1, 62 - 5, 1, 12 - 5, 1,      1, 1,      1, 62 - 5, 1,

               13 - 5, 1, 1,      1, 1,      63 - 5, 1, 13 - 5, 1, 1,      1, 1,
               63 - 5, 1, 13 - 5, 1, 1,      1,      1, 63 - 5, 1, 13 - 5, 1, 1,
               1,      1, 63 - 5, 1, 13 - 5, 1,      1, 1,      1, 63 - 5, 1,
           }},
          {"alt_row_2",
           {
               11 + 5, 1, 1,      1, 1,      61 + 5, 1, 11 + 5, 1, 1,      1, 1,
               61 + 5, 1, 11 + 5, 1, 1,      1,      1, 61 + 5, 1, 11 + 5, 1, 1,
               1,      1, 61 + 5, 1, 11 + 5, 1,      1, 1,      1, 61 + 5, 1,

               12 + 5, 1, 1,      1, 1,      62 + 5, 1, 12 + 5, 1, 1,      1, 1,
               62 + 5, 1, 12 + 5, 1, 1,      1,      1, 62 + 5, 1, 12 + 5, 1, 1,
               1,      1, 62 + 5, 1, 12 + 5, 1,      1, 1,      1, 62 + 5, 1,

               13 + 5, 1, 1,      1, 1,      63 + 5, 1, 13 + 5, 1, 1,      1, 1,
               63 + 5, 1, 13 + 5, 1, 1,      1,      1, 63 + 5, 1, 13 + 5, 1, 1,
               1,      1, 63 + 5, 1, 13 + 5, 1,      1, 1,      1, 63 + 5, 1,
           }},
          {"base_channels",
           {
               11, 21, 31, 41, 51, 61, 71, 11 - 5, 11 + 5,
               11, 21, 31, 41, 51, 61, 71, 11 - 5, 11 + 5,
               11, 21, 31, 41, 51, 61, 71, 11 - 5, 11 + 5,
               11, 21, 31, 41, 51, 61, 71, 11 - 5, 11 + 5,
               11, 21, 31, 41, 51, 61, 71, 11 - 5, 11 + 5,

               12, 22, 32, 42, 52, 62, 72, 12 - 5, 12 + 5,
               12, 22, 32, 42, 52, 62, 72, 12 - 5, 12 + 5,
               12, 22, 32, 42, 52, 62, 72, 12 - 5, 12 + 5,
               12, 22, 32, 42, 52, 62, 72, 12 - 5, 12 + 5,
               12, 22, 32, 42, 52, 62, 72, 12 - 5, 12 + 5,

               13, 23, 33, 43, 53, 63, 73, 13 - 5, 13 + 5,
               13, 23, 33, 43, 53, 63, 73, 13 - 5, 13 + 5,
               13, 23, 33, 43, 53, 63, 73, 13 - 5, 13 + 5,
               13, 23, 33, 43, 53, 63, 73, 13 - 5, 13 + 5,
               13, 23, 33, 43, 53, 63, 73, 13 - 5, 13 + 5,
           }},
          {"diff_channels",
           {
               11, 21, 31, 41, 51, 61, 71, 61 - 5, 61 + 5,
               11, 21, 31, 41, 51, 61, 71, 61 - 5, 61 + 5,
               11, 21, 31, 41, 51, 61, 71, 61 - 5, 61 + 5,
               11, 21, 31, 41, 51, 61, 71, 61 - 5, 61 + 5,
               11, 21, 31, 41, 51, 61, 71, 61 - 5, 61 + 5,

               12, 22, 32, 42, 52, 62, 72, 62 - 5, 62 + 5,
               12, 22, 32, 42, 52, 62, 72, 62 - 5, 62 + 5,
               12, 22, 32, 42, 52, 62, 72, 62 - 5, 62 + 5,
               12, 22, 32, 42, 52, 62, 72, 62 - 5, 62 + 5,
               12, 22, 32, 42, 52, 62, 72, 62 - 5, 62 + 5,

               13, 23, 33, 43, 53, 63, 73, 63 - 5, 63 + 5,
               13, 23, 33, 43, 53, 63, 73, 63 - 5, 63 + 5,
               13, 23, 33, 43, 53, 63, 73, 63 - 5, 63 + 5,
               13, 23, 33, 43, 53, 63, 73, 63 - 5, 63 + 5,
               13, 23, 33, 43, 53, 63, 73, 63 - 5, 63 + 5,
           }}};

  const FillPileupArrayTestData& param = GetParam();
  int size = input.size() * param.expected_num_rows * width *
             param.expected_num_channels;
  std::vector<unsigned char> pileup_image(size, 0);
  FillPileupArray(input, alt_rows, param.alt_aligned_pileup, &pileup_image,
                  size, 0, param.alt_image_row_indices);
  std::vector<unsigned char> expected_pilup_image = ConcatExpectedPileupImages(
      param.expected_pileup_images, expected_pileup_image_by_name);
  EXPECT_THAT(pileup_image, testing::ElementsAreArray(expected_pilup_image));
}

INSTANTIATE_TEST_SUITE_P(
    FillPileupArrayTests, FillPileupArrayTest,
    testing::ValuesIn(std::vector<FillPileupArrayTestData>({
        {
            .test_name = "test_alt_aligned_pileup_none",
            .expected_num_rows = 1,
            .expected_num_channels = 7,
            .alt_aligned_pileup = AltAlignedPileup::kNone,
            .alt_image_row_indices = {},
            .expected_pileup_images = {"ref_row"},
        },
        {.test_name = "test_alt_aligned_pileup_base_channels",
         .expected_num_rows = 1,
         .expected_num_channels = 9,
         .alt_aligned_pileup = AltAlignedPileup::kBaseChannels,
         .alt_image_row_indices = {},
         .expected_pileup_images = {"base_channels"}},
        {
            .test_name = "test_alt_aligned_pileup_diff_channels",
            .expected_num_rows = 1,
            .expected_num_channels = 9,
            .alt_aligned_pileup = AltAlignedPileup::kDiffChannels,
            .alt_image_row_indices = {},
            .expected_pileup_images = {"diff_channels"},
        },
        {.test_name = "test_alt_aligned_pileup_rows",
         .expected_num_rows = 3,
         .expected_num_channels = 7,
         .alt_aligned_pileup = AltAlignedPileup::kRows,
         .alt_image_row_indices = {0, 1},
         .expected_pileup_images = {"ref_row", "alt_row_1", "alt_row_2"}},
        {
            .test_name = "test_alt_aligned_pileup_single_row_first",
            .expected_num_rows = 2,
            .expected_num_channels = 7,
            .alt_aligned_pileup = AltAlignedPileup::kSingleRow,
            .alt_image_row_indices = {0},
            .expected_pileup_images = {"ref_row", "alt_row_1"},
        },
        {
            .test_name = "test_alt_aligned_pileup_single_row_second",
            .expected_num_rows = 2,
            .expected_num_channels = 7,
            .alt_aligned_pileup = AltAlignedPileup::kSingleRow,
            .alt_image_row_indices = {1},
            .expected_pileup_images = {"ref_row", "alt_row_2"},
        },
    })));

struct FillPileupArrayBySampleTestData {
  std::string test_name;
  std::string alt_aligned_pileup;
  std::string sample_1_alt_aligned_pileup;
  std::string sample_2_alt_aligned_pileup;
  std::vector<std::string> alt_combination;
  std::vector<std::string> expected_pileup_images;
};

class FillPileupArrayBySampleTest
    : public testing::TestWithParam<FillPileupArrayBySampleTestData> {};

INSTANTIATE_TEST_SUITE_P(
    FillPileupArrayBySampleTests, FillPileupArrayBySampleTest,
    testing::ValuesIn(std::vector<FillPileupArrayBySampleTestData>(
        {{.test_name = "test_alt_aligned_pileup_global_none",
          .alt_aligned_pileup = "none",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_2_ref_row"}},
         {.test_name = "test_alt_aligned_pileup_global_rows",
          .alt_aligned_pileup = "rows",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_1_alt_row_1",
                                     "sample_1_alt_row_2", "sample_2_ref_row",
                                     "sample_2_alt_row_1",
                                     "sample_2_alt_row_2"}},
         {.test_name = "test_alt_aligned_pileup_global_single_row",
          .alt_aligned_pileup = "single_row",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_1_alt_row_1",
                                     "sample_2_ref_row", "sample_2_alt_row_1"}},
         {.test_name =
              "test_alt_aligned_pileup_global_single_row_second_is_longer",
          .alt_aligned_pileup = "single_row",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A", "AA"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_1_alt_row_2",
                                     "sample_2_ref_row", "sample_2_alt_row_2"}},
         {.test_name = "test_alt_aligned_pileup_global_rows_sample_1_none",
          .alt_aligned_pileup = "rows",
          .sample_1_alt_aligned_pileup = "none",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_2_ref_row",
                                     "sample_2_alt_row_1",
                                     "sample_2_alt_row_2"}},
         {.test_name = "test_alt_aligned_pileup_global_rows_sample_2_none",
          .alt_aligned_pileup = "rows",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "none",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_1_alt_row_1",
                                     "sample_1_alt_row_2", "sample_2_ref_row"}},
         {.test_name =
              "test_alt_aligned_pileup_global_single_row_sample_1_none",
          .alt_aligned_pileup = "single_row",
          .sample_1_alt_aligned_pileup = "none",
          .sample_2_alt_aligned_pileup = "",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_2_ref_row",
                                     "sample_2_alt_row_1"}},
         {.test_name =
              "test_alt_aligned_pileup_global_single_row_sample_2_none",
          .alt_aligned_pileup = "single_row",
          .sample_1_alt_aligned_pileup = "",
          .sample_2_alt_aligned_pileup = "none",
          .alt_combination = {"A"},
          .expected_pileup_images = {"sample_1_ref_row", "sample_1_alt_row_1",
                                     "sample_2_ref_row"}}})));

TEST_P(FillPileupArrayBySampleTest, FillPileupArrayBySampleTestCases) {
  int height = 2;
  int width = 5;
  int num_channels = 2;

  std::vector<std::unique_ptr<ImageRow>> sample_1_input;
  sample_1_input.emplace_back(MakeImageRow(
      {{11, 11, 11, 11, 11}, {21, 21, 21, 21, 21}}, width, num_channels));
  sample_1_input.emplace_back(MakeImageRow(
      {{12, 12, 12, 12, 12}, {22, 22, 22, 22, 22}}, width, num_channels));
  // In order to distinguish channel 8 and 9 "-5" is added to channel 8 and
  // "+5" is added to channel 9.
  std::vector<std::vector<std::unique_ptr<ImageRow>>> sample_1_alt_rows(2);
  sample_1_alt_rows[0].emplace_back(
      MakeImageRow({{11 - 5, 11 - 5, 11 - 5, 11 - 5, 11 - 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_1_alt_rows[0].emplace_back(
      MakeImageRow({{12 - 5, 12 - 5, 12 - 5, 12 - 5, 12 - 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_1_alt_rows[1].emplace_back(
      MakeImageRow({{11 + 5, 11 + 5, 11 + 5, 11 + 5, 11 + 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_1_alt_rows[1].emplace_back(
      MakeImageRow({{12 + 5, 12 + 5, 12 + 5, 12 + 5, 12 + 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));

  std::vector<std::unique_ptr<ImageRow>> sample_2_input;
  sample_2_input.emplace_back(MakeImageRow(
      {{51, 51, 51, 51, 51}, {61, 61, 61, 61, 61}}, width, num_channels));
  sample_2_input.emplace_back(MakeImageRow(
      {{52, 52, 52, 52, 52}, {62, 62, 62, 62, 62}}, width, num_channels));
  // In order to distinguish channel 8 and 9 "-5" is added to channel 8 and
  // "+5" is added to channel 9.
  std::vector<std::vector<std::unique_ptr<ImageRow>>> sample_2_alt_rows(2);
  sample_2_alt_rows[0].emplace_back(
      MakeImageRow({{51 - 5, 51 - 5, 51 - 5, 51 - 5, 51 - 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_2_alt_rows[0].emplace_back(
      MakeImageRow({{52 - 5, 52 - 5, 52 - 5, 52 - 5, 52 - 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_2_alt_rows[1].emplace_back(
      MakeImageRow({{51 + 5, 51 + 5, 51 + 5, 51 + 5, 51 + 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));
  sample_2_alt_rows[1].emplace_back(
      MakeImageRow({{52 + 5, 52 + 5, 52 + 5, 52 + 5, 52 + 5}, {1, 1, 1, 1, 1}},
                   width, num_channels));

  std::vector<std::vector<std::unique_ptr<ImageRow>>> image_per_sample;
  image_per_sample.push_back(std::move(sample_1_input));
  image_per_sample.push_back(std::move(sample_2_input));

  std::vector<std::vector<std::vector<std::unique_ptr<ImageRow>>>>
      alt_image_per_sample;
  alt_image_per_sample.push_back(std::move(sample_1_alt_rows));
  alt_image_per_sample.push_back(std::move(sample_2_alt_rows));

  std::unordered_map<std::string, std::vector<unsigned char>>
      expected_pileup_image_by_name = {
          {"sample_1_ref_row",
           {
               11, 21, 11, 21, 11, 21, 11, 21, 11, 21,
               12, 22, 12, 22, 12, 22, 12, 22, 12, 22,
           }},
          {"sample_1_alt_row_1",
           {
               11 - 5, 1, 11 - 5, 1, 11 - 5, 1, 11 - 5, 1, 11 - 5, 1,
               12 - 5, 1, 12 - 5, 1, 12 - 5, 1, 12 - 5, 1, 12 - 5, 1,
           }},
          {"sample_1_alt_row_2",
           {
               11 + 5, 1, 11 + 5, 1, 11 + 5, 1, 11 + 5, 1, 11 + 5, 1,
               12 + 5, 1, 12 + 5, 1, 12 + 5, 1, 12 + 5, 1, 12 + 5, 1,
           }},
          {"sample_2_ref_row",
           {
               51, 61, 51, 61, 51, 61, 51, 61, 51, 61,
               52, 62, 52, 62, 52, 62, 52, 62, 52, 62,
           }},
          {"sample_2_alt_row_1",
           {
               51 - 5, 1, 51 - 5, 1, 51 - 5, 1, 51 - 5, 1, 51 - 5, 1,

               52 - 5, 1, 52 - 5, 1, 52 - 5, 1, 52 - 5, 1, 52 - 5, 1,
           }},
          {"sample_2_alt_row_2",
           {
               51 + 5, 1, 51 + 5, 1, 51 + 5, 1, 51 + 5, 1, 51 + 5, 1,

               52 + 5, 1, 52 + 5, 1, 52 + 5, 1, 52 + 5, 1, 52 + 5, 1,
           }}};

  const FillPileupArrayBySampleTestData& param = GetParam();
  MakeExamplesOptions options;
  options.mutable_pic_options()->set_height(height);
  options.mutable_pic_options()->set_width(width);
  options.mutable_pic_options()->set_num_channels(num_channels);
  options.mutable_pic_options()->set_alt_aligned_pileup(
      param.alt_aligned_pileup);
  options.mutable_pic_options()->set_types_to_alt_align("all");

  // Add first sample
  SampleOptions* sample_options = options.add_sample_options();
  sample_options->set_name("sample_1");
  sample_options->set_pileup_height(height);
  sample_options->set_alt_aligned_pileup(param.sample_1_alt_aligned_pileup);

  // Add second sample
  SampleOptions* sample_options_2 = options.add_sample_options();
  sample_options_2->set_name("sample_2");
  sample_options_2->set_pileup_height(height);
  sample_options_2->set_alt_aligned_pileup(param.sample_2_alt_aligned_pileup);

  int pileup_image_height = CalculatePileupImageHeight(options);
  int size = pileup_image_height * width * num_channels;
  std::vector<unsigned char> pileup_image(size, 0);
  FillPileupArrayBySample(image_per_sample, alt_image_per_sample, options,
                          param.alt_combination, &pileup_image, size);
  std::vector<unsigned char> expected_pilup_image = ConcatExpectedPileupImages(
      param.expected_pileup_images, expected_pileup_image_by_name);
  EXPECT_EQ(pileup_image.size(), size);
  EXPECT_THAT(pileup_image, testing::ElementsAreArray(expected_pilup_image));
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
