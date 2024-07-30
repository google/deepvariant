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

#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include "deepvariant/variant_calling_multisample.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/container/node_hash_map.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace multi_sample {

using absl::StrCat;
using nucleus::EqualsProto;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

constexpr char kSampleName[] = "MySampleName";
constexpr char kChr[] = "chr1";
constexpr int64_t kStart = 10;

AlleleCount MakeTestAlleleCount(int total_n, int alt_n,
                                const std::string& sample_id,
                                const std::string& ref = "A",
                                const std::string& alt = "C", int start = 100) {
  CHECK_GE(total_n, alt_n) << "Total number of reads must be >= n alt reads";
  AlleleCount allele_count;
  *(allele_count.mutable_position()) = nucleus::MakePosition("chr1", start);
  allele_count.set_ref_base(ref);
  allele_count.set_ref_supporting_read_count(total_n - alt_n);
  const Allele read_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, 1);
  for (int i = 0; i < alt_n; ++i) {
    (*allele_count
          .mutable_read_alleles())[StrCat(sample_id, "_read_", i)] =
        read_allele;

    Allele* new_allele =
        (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
    *new_allele = read_allele;
  }
  return allele_count;
}

VariantCallerOptions BasicOptions() {
  // Set basic options to avoid premature test failures.
  VariantCallerOptions options;
  options.set_sample_name(kSampleName);
  options.set_ploidy(2);

  return options;
}

Variant MakeExpectedVariant(const std::string& ref,
                            const std::vector<std::string>& alts,
                            const int64_t start = kStart) {
  Variant variant;
  variant.set_reference_name(kChr);
  variant.set_start(start);
  variant.set_reference_bases(ref);
  for (const std::string& alt_allele : alts)
    variant.add_alternate_bases(alt_allele);

  if (alts.empty()) {
    variant.set_end(variant.start() + 1);
    variant.add_alternate_bases(kGVCFAltAllele);
    CHECK(google::protobuf::TextFormat::ParseFromString(
        "genotype: 0 genotype: 0 "
        "genotype_likelihood: -0.47712125472 "
        "genotype_likelihood: -0.47712125472 "
        "genotype_likelihood: -0.47712125472 "
        "info: { key: \"GQ\" value { values { int_value: 1 } } }",
        variant.add_calls()));
  } else {
    // End is start + ref length according to Variant.proto spec.
    variant.set_end(variant.start() + ref.length());
    CHECK(google::protobuf::TextFormat::ParseFromString("genotype: -1 genotype: -1",
                                              variant.add_calls()));
  }

  VariantCall* call = variant.mutable_calls(0);
  call->set_call_set_name(kSampleName);
  return variant;
}

Variant WithCounts(const Variant& base_variant, const std::vector<int>& ad,
                   int dp = -1) {
  CHECK(!ad.empty() || dp != -1) << "Either AD or DP must be provided.";
  Variant variant(base_variant);
  VariantCall* call = variant.mutable_calls(0);
  if (ad.empty()) {
    nucleus::SetInfoField(kDPFormatField, dp, call);
  } else {
    if (dp == -1) dp = std::accumulate(ad.begin(), ad.end(), 0);
    nucleus::SetInfoField(kDPFormatField, dp, call);
    nucleus::SetInfoField(kADFormatField, ad, call);
    std::vector<double> vaf;
    // Skip the first one in ad which is ref.
    for (size_t i = 1; i < ad.size(); ++i) {
      vaf.push_back(1.0 * ad[i] / dp);
    }
    nucleus::SetInfoField(kVAFFormatField, vaf, call);
  }
  return variant;
}

// Test max_fraction_snps_for_non_target_sample.
TEST(VariantCallingSomaticTest, TestCallVariantWithMaxFractionForNormal) {
  // Tumor:
  // SNP A -> T with 19 reads support and 1 reads ref support
  //
  // Normal:
  // SNP A -> T with 7 reads support and 3 reads ref support
  //
  // set_min_fraction_snps >= 0.1

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["tumor"] = MakeTestAlleleCount(20, 19, "tumor", "A", "T", 10);
  allele_counts["normal"] = MakeTestAlleleCount(10, 7, "normal", "A", "T", 10);

  VariantCallerOptions options = BasicOptions();
  options.set_min_fraction_snps(0.1);
  options.set_max_fraction_snps_for_non_target_sample(0.0);  // 0 means not set.
  // Set min_fraction_multiplier because we want to test the `CallVariant` for
  // somatic mode.
  options.set_min_fraction_multiplier(std::numeric_limits<float>::infinity());
  const VariantCaller caller(options);
  const std::optional<DeepVariantCall> optional_variant =
      caller.CallVariant(allele_counts, "tumor");
  EXPECT_TRUE(static_cast<bool>(optional_variant));
  Variant variant = WithCounts(MakeExpectedVariant("A", {"T"}, 10), {1, 19});
  EXPECT_THAT(optional_variant->variant(), EqualsProto(variant));

  const double EPSILON = 1e-6;

  // Now we set max_fraction_snps_for_non_target_sample to 0.7+EPSILON, which
  // should still create the variant.
  options.set_max_fraction_snps_for_non_target_sample(0.7+EPSILON);
  const VariantCaller caller2(options);
  const std::optional<DeepVariantCall> optional_variant2 =
      caller2.CallVariant(allele_counts, "tumor");
  EXPECT_TRUE(static_cast<bool>(optional_variant2));
  EXPECT_THAT(optional_variant2->variant(), EqualsProto(variant));

  // Now we set max_fraction_snps_for_non_target_sample to 0.7-EPSILON, which
  // should stop the variant from being created.
  options.set_max_fraction_snps_for_non_target_sample(0.7-EPSILON);
  const VariantCaller caller3(options);
  const std::optional<DeepVariantCall> optional_variant3 =
      caller3.CallVariant(allele_counts, "tumor");
  EXPECT_FALSE(static_cast<bool>(optional_variant3));
}

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
