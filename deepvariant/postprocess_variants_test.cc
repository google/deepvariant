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

#include "deepvariant/postprocess_variants.h"

#include <fstream>
#include <string>
#include <vector>

// #include "google/protobuf/struct.pb.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/io/compression.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::Variant;
using ::testing::ElementsAreArray;

namespace {

CallVariantsOutput CreateSingleSiteCalls(absl::string_view reference_name,
                                         int start,
                                         int end) {
  CallVariantsOutput single_site_call;
  // Add one call to fulfill the assumption of variant having one call.
  single_site_call.mutable_variant()->add_calls();
  single_site_call.mutable_variant()->set_reference_name(
      string(reference_name));
  single_site_call.mutable_variant()->set_start(start);
  single_site_call.mutable_variant()->set_end(end);
  return single_site_call;
}

CallVariantsOutput CreateSingleSiteCalls(absl::string_view reference_name,
                                         int start,
                                         int end, double quality) {
  CallVariantsOutput single_site_call =
      CreateSingleSiteCalls(reference_name, start, end);
  single_site_call.mutable_variant()->set_quality(quality);
  return single_site_call;
}

// region Phasing tests helpers
void WriteSwitchesFile(const std::string& path,
                       absl::Span<const std::string> lines) {
  std::ofstream outfile(path);
  for (const auto& line : lines) {
    outfile << line << "\n";
  }
  outfile.close();
}

Variant CreatePhasedVariant(std::string chrom, int start, int end,
                            std::string ps_contig,
                            bool is_first_variant_in_phase_set, bool is_phased,
                            absl::Span<const int> gt, bool add_alt_ps = true) {
  Variant variant;
  variant.set_reference_name(chrom);
  variant.set_start(start);
  variant.set_end(end);
  variant.add_alternate_bases("T");
  variant.set_reference_bases("C");
  auto* call = variant.add_calls();
  call->set_is_phased(is_phased);
  for (int g : gt) {
    call->add_genotype(g);
  }

  if (!ps_contig.empty()) {
    (*variant.mutable_info())[kPhaseSetContigTag]
        .add_values()
        ->set_string_value(ps_contig);
  }
  if (is_first_variant_in_phase_set) {
    (*variant.mutable_info())[kFirstVariantInPhaseSetTag]
        .add_values()
        ->set_bool_value(true);
  }
  if (add_alt_ps) {
    auto& list_value = (*variant.mutable_info())[kAltPhaseSetTag];
    list_value.add_values()->set_number_value(1);  // dummy ALT_PS
    list_value.add_values()->set_number_value(2);
  }
  return variant;
}

// endregion

}  // namespace

TEST(ProcessSingleSiteCallTfRecords, BasicCase) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 1000, 1001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.9));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.7));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path,
                                 std::vector<nucleus::genomics::v1::Range>());
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 5);
  EXPECT_EQ(output[0].variant().reference_name(), "chr1");
  EXPECT_EQ(output[1].variant().reference_name(), "chr10");
  EXPECT_EQ(output[2].variant().reference_name(), "chr10");
  EXPECT_EQ(output[3].variant().reference_name(), "chr10");
  EXPECT_EQ(output[4].variant().reference_name(), "chr10");
  EXPECT_EQ(output[0].variant().start(), 1);
  EXPECT_EQ(output[1].variant().start(), 1000);
  EXPECT_EQ(output[2].variant().start(), 2000);
  EXPECT_EQ(output[3].variant().start(), 2000);
  EXPECT_EQ(output[4].variant().start(), 2000);
  EXPECT_EQ(output[0].variant().end(), 2);
  EXPECT_EQ(output[1].variant().end(), 1001);
  EXPECT_EQ(output[2].variant().end(), 2001);
  EXPECT_EQ(output[3].variant().end(), 2002);
  EXPECT_EQ(output[4].variant().end(), 2002);
  // Order of calls with the same reference, start, end should be preserved.
  EXPECT_EQ(output[3].variant().quality(), 0.9);
  EXPECT_EQ(output[4].variant().quality(), 0.7);
}

TEST(ProcessSingleSiteCallTfRecords, SplitByRange) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  std::vector<nucleus::genomics::v1::Range> ranges = {
      nucleus::MakeRange("chr1", 0, 100)};
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 1000, 1001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.9));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.7));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path, ranges);
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 1);
  EXPECT_EQ(output[0].variant().reference_name(), "chr1");
  EXPECT_EQ(output[0].variant().start(), 1);
  EXPECT_EQ(output[0].variant().end(), 2);
}

TEST(ProcessSingleSiteCallTfRecords, SplitByRangesHandleNeighbors) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  std::vector<nucleus::genomics::v1::Range> ranges = {
      nucleus::MakeRange("chr10", 0, 100),
      nucleus::MakeRange("chr10", 100, 200),
  };
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 50, 51));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 100, 110));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 199, 200));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 199, 210));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 200, 210));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path, ranges);
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 4);
  EXPECT_EQ(output[0].variant().reference_name(), "chr10");
  EXPECT_EQ(output[0].variant().start(), 50);
  EXPECT_EQ(output[0].variant().end(), 51);
  EXPECT_EQ(output[1].variant().reference_name(), "chr10");
  EXPECT_EQ(output[1].variant().start(), 100);
  EXPECT_EQ(output[1].variant().end(), 110);
  EXPECT_EQ(output[2].variant().reference_name(), "chr10");
  EXPECT_EQ(output[2].variant().start(), 199);
  EXPECT_EQ(output[2].variant().end(), 200);
  EXPECT_EQ(output[3].variant().reference_name(), "chr10");
  EXPECT_EQ(output[3].variant().start(), 199);
  EXPECT_EQ(output[3].variant().end(), 210);
}

// region Phasing tests

struct MaybeSwapPhaseTestData {
  std::string name;
  VariantPhaseInformation phase_info;
  bool is_phased;
  std::vector<int> initial_gt;
  std::vector<int> expected_gt;
  bool expect_ps_tag;
};

class MaybeSwapPhaseTest
    : public ::testing::TestWithParam<MaybeSwapPhaseTestData> {};

TEST_P(MaybeSwapPhaseTest, MaybeSwapPhaseTestCases) {
  const MaybeSwapPhaseTestData& param = GetParam();
  Variant variant = CreatePhasedVariant("chr1", 110, 111, "1_1", false,
                                        param.is_phased, param.initial_gt);

  MaybeSwapPhase(&variant, param.phase_info);

  EXPECT_THAT(variant.calls(0).genotype(), ElementsAreArray(param.expected_gt));
  if (param.expect_ps_tag) {
    ASSERT_TRUE(variant.calls(0).info().contains(kPhaseSetTag));
    EXPECT_EQ(variant.calls(0).info().at(kPhaseSetTag).values(0).int_value(),
              param.phase_info.first_variant_in_phase_set_start_position + 1);
  } else {
    EXPECT_FALSE(variant.calls(0).info().contains(kPhaseSetTag));
  }
}

INSTANTIATE_TEST_SUITE_P(
    MaybeSwapPhaseTests, MaybeSwapPhaseTest,
    ::testing::ValuesIn(std::vector<MaybeSwapPhaseTestData>({
        {"NoSwapOnMatch",
         {"1", "1", PhaseSetStitchingStatus::MATCH, false, 100, false},
         true,
         {0, 1},
         {0, 1},
         true},
        {"SwapOnSwitch",
         {"1", "1", PhaseSetStitchingStatus::SWITCH, false, 100, false},
         true,
         {0, 1},
         {1, 0},
         true},
        {"NoSwapOnNotEnoughOverlap",
         {"1", "1", PhaseSetStitchingStatus::NOT_ENOUGH_OVERLAP, false, 100,
          false},
         true,
         {0, 1},
         {0, 1},
         true},
        {"NoSwapIfHom",
         {"1", "1", PhaseSetStitchingStatus::SWITCH, false, 100, false},
         true,
         {1, 1},
         {1, 1},
         true},
        {"NoActionIfNotPhased",
         {"1", "1", PhaseSetStitchingStatus::SWITCH, false, 100, false},
         false,
         {0, 1},
         {0, 1},
         false},
    })),
    [](const ::testing::TestParamInfo<MaybeSwapPhaseTest::ParamType>& info) {
      return info.param.name;
    });

struct StitchPhaseSetsTestData {
  std::string name;
  std::vector<Variant> input_variants;
  std::vector<std::vector<int>> expected_gts;
  std::vector<int> expected_ps_tags;
};

class StitchPhaseSetsTest
    : public ::testing::TestWithParam<StitchPhaseSetsTestData> {};

TEST_P(StitchPhaseSetsTest, StitchPhaseSetsTestCases) {
  const StitchPhaseSetsTestData& param = GetParam();
  const std::string switches_path = nucleus::MakeTempFile("switches.tsv");
  /* shard 0: MATCH, shard 1: SWITCH, shard 2: NOT_ENOUGH_OVERLAP */
  WriteSwitchesFile(switches_path,
                    {"0\t0\t0", "0\t1\t1", "0\t2\t2", "1\t0\t0"});
  const std::string input_tfrecord_path =
      nucleus::MakeTempFile("stitch.in.tfrecord");
  const std::string output_tfrecord_path =
      nucleus::MakeTempFile("stitch.out.tfrecord");
  nucleus::WriteProtosToTFRecord(param.input_variants, input_tfrecord_path,
                                 tensorflow::io::compression::kNone);

  StitchPhaseSets({input_tfrecord_path}, switches_path, {output_tfrecord_path});

  std::vector<Variant> output =
      nucleus::ReadProtosFromTFRecord<Variant>(output_tfrecord_path);
  ASSERT_EQ(output.size(), param.input_variants.size());
  for (int i = 0; i < output.size(); ++i) {
    EXPECT_THAT(output[i].calls(0).genotype(),
                ElementsAreArray(param.expected_gts[i]));
    EXPECT_EQ(output[i].calls(0).info().at(kPhaseSetTag).values(0).int_value(),
              param.expected_ps_tags[i]);
  }
}

INSTANTIATE_TEST_SUITE_P(
    StitchPhaseSetsTests, StitchPhaseSetsTest,
    ::testing::ValuesIn(std::vector<StitchPhaseSetsTestData>({
        {
            "StichAllVariantsIntoOnePhaseSet",
            {
                CreatePhasedVariant("chr1", 100, 101, "0-0", true, true,
                                    {0, 1}),
                CreatePhasedVariant("chr1", 110, 111, "0-0", false, true,
                                    {0, 1}),
                CreatePhasedVariant("chr1", 200, 201, "0-1", false, true,
                                    {0, 1}),
                CreatePhasedVariant("chr1", 210, 211, "0-1", false, true,
                                    {0, 1}),
            },
            {{0, 1}, {0, 1}, {1, 0}, {1, 0}},
            {101, 101, 101, 101, 101, 101},
        },
        {"NotEnoughOverlapStartsNewPhaseSet",
         {
             CreatePhasedVariant("chr1", 100, 101, "0-0", true, true, {0, 1}),
             CreatePhasedVariant("chr1", 110, 111, "0-0", false, true, {0, 1}),
             CreatePhasedVariant("chr1", 200, 201, "0-2", false, true, {0, 1}),
             CreatePhasedVariant("chr1", 210, 211, "0-2", false, true, {0, 1}),
         },
         {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
         {101, 101, 201, 201}},
        {"FirstVariantInPhaseSetStartsNewPhaseSet",
         {
             CreatePhasedVariant("chr1", 100, 101, "0-0", true, true, {0, 1}),
             CreatePhasedVariant("chr1", 110, 111, "0-0", false, true, {0, 1}),
             CreatePhasedVariant("chr1", 200, 201, "0-1", true, true, {0, 1}),
             CreatePhasedVariant("chr1", 210, 211, "0-1", false, true, {0, 1}),
         },
         {{0, 1}, {0, 1}, {1, 0}, {1, 0}},
         {101, 101, 201, 201}},
        {"NewRegionWithFirstVariantInPhaseSetStartsNewPhaseSet",
         {
             CreatePhasedVariant("chr1", 100, 101, "0-0", true, true, {0, 1}),
             CreatePhasedVariant("chr1", 110, 111, "0-0", false, true, {0, 1}),
             CreatePhasedVariant("chr1", 200, 201, "1-0", true, true, {0, 1}),
             CreatePhasedVariant("chr1", 210, 211, "1-0", false, true, {0, 1}),
         },
         {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
         {101, 101, 201, 201}},
    })),
    [](const ::testing::TestParamInfo<StitchPhaseSetsTest::ParamType>& info) {
      return info.param.name;
    });

// endregion

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
