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

#include "deepvariant/variant_calling_multisample.h"

#include <optional>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace multi_sample {

using absl::StrCat;

constexpr char kSampleName[] = "MySampleName";

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
    (*allele_count.mutable_read_alleles())[StrCat(sample_id, "_read_", i)] =
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

// Test small_model_vaf_context_window_size = 5.
TEST(VariantCallingTest,
     TestCallVariantAddAdjacentAlleleFractionsAtPositionSize5) {
  AlleleCount allele_count =
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 10);
  const std::vector<AlleleCount>& target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "A", "T", 7),
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 8),
      MakeTestAlleleCount(20, 17, "sample", "A", "T", 9),
      allele_count,
      MakeTestAlleleCount(20, 20, "sample", "A", "T", 11),
      MakeTestAlleleCount(20, 0, "sample", "A", "T", 12),
      MakeTestAlleleCount(20, 10, "sample", "A", "T", 13),
  };
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;
  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] = target_sample_allele_counts.begin() + 3;

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  int window_size = 5;
  VariantCallerOptions options = BasicOptions();
  options.set_small_model_vaf_context_window_size(window_size);

  const VariantCaller caller(options);
  const std::optional<DeepVariantCall> optional_variant =
      caller.CallVariant(allele_counts, "sample", &target_sample_allele_counts,
                         &allele_counter_iterators["sample"]);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
  EXPECT_THAT(call.allele_frequency_at_position().at(8), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(9), 85);
  EXPECT_THAT(call.allele_frequency_at_position().at(10), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(11), 100);
  EXPECT_THAT(call.allele_frequency_at_position().at(12), 0);
}

// Test small_model_vaf_context_window_size = 3.
TEST(VariantCallingTest,
     TestCallVariantAddAdjacentAlleleFractionsAtPositionSize3) {
  AlleleCount allele_count =
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 10);
  const std::vector<AlleleCount>& target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "A", "T", 7),
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 8),
      MakeTestAlleleCount(20, 17, "sample", "A", "T", 9),
      allele_count,
      MakeTestAlleleCount(20, 20, "sample", "A", "T", 11),
      MakeTestAlleleCount(20, 0, "sample", "A", "T", 12),
      MakeTestAlleleCount(20, 10, "sample", "A", "T", 13),
  };
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;
  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] = target_sample_allele_counts.begin() + 3;

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  int window_size = 3;
  VariantCallerOptions options = BasicOptions();
  options.set_small_model_vaf_context_window_size(window_size);

  const VariantCaller caller(options);
  const std::optional<DeepVariantCall> optional_variant =
      caller.CallVariant(allele_counts, "sample", &target_sample_allele_counts,
                         &allele_counter_iterators["sample"]);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
  EXPECT_THAT(call.allele_frequency_at_position().at(9), 85);
  EXPECT_THAT(call.allele_frequency_at_position().at(10), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(11), 100);
}

// Test small_model_vaf_context_window_size = 0.
TEST(VariantCallingTest,
     TestCallVariantAddAdjacentAlleleFractionsAtPositionSize0) {
  AlleleCount allele_count =
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 10);
  const std::vector<AlleleCount>& target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "A", "T", 7),
      MakeTestAlleleCount(20, 19, "sample", "A", "T", 8),
      MakeTestAlleleCount(20, 17, "sample", "A", "T", 9),
      allele_count,
      MakeTestAlleleCount(20, 20, "sample", "A", "T", 11),
      MakeTestAlleleCount(20, 0, "sample", "A", "T", 12),
      MakeTestAlleleCount(20, 10, "sample", "A", "T", 13),
  };
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;
  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] = target_sample_allele_counts.begin() + 3;

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  int window_size = 0;
  VariantCallerOptions options = BasicOptions();
  options.set_small_model_vaf_context_window_size(window_size);

  const VariantCaller caller(options);
  const std::optional<DeepVariantCall> optional_variant =
      caller.CallVariant(allele_counts, "sample", &target_sample_allele_counts,
                         &allele_counter_iterators["sample"]);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
}

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
