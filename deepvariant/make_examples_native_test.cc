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

#include "deepvariant/make_examples_native.h"

#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/log/check.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "google/protobuf/text_format.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Variant = nucleus::genomics::v1::Variant;
using VariantCall = nucleus::genomics::v1::VariantCall;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

constexpr char kSampleName[] = "MySampleName";
constexpr char kChr[] = "chr1";
constexpr int64_t kStart = 10;

Variant MakeVariant(absl::string_view ref,
                    const std::vector<absl::string_view>& alts,
                    const int64_t start = kStart) {
  Variant variant;
  variant.set_reference_name(kChr);
  variant.set_start(start);
  variant.set_reference_bases(ref.data(), ref.size());
  for (const auto alt_allele : alts)
    variant.add_alternate_bases(alt_allele.data(), alt_allele.size());

  // End is start + ref length according to Variant.proto spec.
  variant.set_end(variant.start() + ref.length());
  CHECK(google::protobuf::TextFormat::ParseFromString("genotype: -1 genotype: -1",
                                            variant.add_calls()));

  VariantCall* call = variant.mutable_calls(0);
  call->set_call_set_name(kSampleName);
  return variant;
}

struct AltAlleleCombinationsTestData {
  PileupImageOptions_MultiAllelicMode mode;
  Variant variant;
  std::vector<std::vector<std::string>> expected_allele_combinations;
};

class AltAlleleCombinationsTest
    : public testing::TestWithParam<AltAlleleCombinationsTestData> {};

TEST_P(AltAlleleCombinationsTest, AltAlleleCombinationsTestCases) {
  const AltAlleleCombinationsTestData& param = GetParam();
  MakeExamplesOptions options;
  options.mutable_pic_options()->set_multi_allelic_mode(param.mode);
  ExamplesGenerator generator(options, /*test_mode=*/true);

  EXPECT_THAT(
      ExamplesGeneratorPeer::CallAltAlleleCombinations(generator,
                                                       param.variant),
      UnorderedElementsAreArray(param.expected_allele_combinations));
}

INSTANTIATE_TEST_SUITE_P(
    AltAlleleCombinationsTests, AltAlleleCombinationsTest,
    ValuesIn(std::vector<AltAlleleCombinationsTestData>({
        {.mode = deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES,
         .variant = MakeVariant("A", {"T"}),
         .expected_allele_combinations = {{"T"}}},
        {.mode = deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES,
         .variant = MakeVariant("AT", {"A"}),
         .expected_allele_combinations = {{"A"}}},
        {.mode = deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES,
         .variant = MakeVariant("A", {"ATT"}),
         .expected_allele_combinations = {{"ATT"}}},
        {.mode = deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES,
        .variant = MakeVariant("AT", {"A", "ATT"}),
         .expected_allele_combinations = {{"A"}, {"ATT"}, {"A", "ATT"}}},
        {.mode = deepvariant::PileupImageOptions::NO_HET_ALT_IMAGES,
        .variant = MakeVariant("AT", {"A", "ATT"}),
         .expected_allele_combinations = {{"A"}, {"ATT"}}},
    })));

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
