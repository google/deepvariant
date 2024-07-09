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

#include <memory>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/testing_utils.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/log/log.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Read = nucleus::genomics::v1::Read;

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

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
