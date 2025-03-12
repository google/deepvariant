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

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/testing_utils.h"
#include "deepvariant/utils.h"
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace multi_sample {

using absl::StrCat;
using ContigInfo = nucleus::genomics::v1::ContigInfo;
using ReferenceSequence = nucleus::genomics::v1::ReferenceSequence;

constexpr char kSampleName[] = "MySampleName";

AlleleCount MakeTestAlleleCount(
    int total_n, int alt_n, const std::string& sample_id,
    const std::string& ref = "A", const std::string& alt = "C", int start = 100,
    int methylated_ref_n = 0, int methylated_alt_n = 0,
    bool track_ref_reads = false) {
  CHECK_GE(total_n, alt_n) << "Total number of reads must be ≥ n alt reads";
  CHECK_GE(total_n - alt_n, methylated_ref_n)
      << "Methylated ref reads must be ≤ total ref reads";
  CHECK_GE(alt_n, methylated_alt_n)
      << "Methylated alt reads must be ≤ alt reads";

  AlleleCount allele_count;
  *(allele_count.mutable_position()) = nucleus::MakePosition("chr1", start);
  allele_count.set_ref_base(ref);
  allele_count.set_ref_supporting_read_count(total_n - alt_n);

  // Create reference alleles if we're tracking ref reads,
  // marking some as methylated
  if (track_ref_reads) {
    for (int i = 0; i < total_n - alt_n; ++i) {
      bool is_methylated = (i < methylated_ref_n);
      Allele ref_allele = MakeAllele(ref, AlleleType::REFERENCE, 1,
                                    false, 30, 30, false, is_methylated);
      (*allele_count.mutable_read_alleles())[StrCat(sample_id,
                                                    "_ref_read_", i)] =
                                                    ref_allele;

      Allele* new_allele =
          (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
      *new_allele = ref_allele;
    }
  }

  // Create alternate alleles, marking some as methylated
  for (int i = 0; i < alt_n; ++i) {
    bool is_methylated = (i < methylated_alt_n);
    Allele alt_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, 1,
                                   false, 30, 30, false, is_methylated);
    (*allele_count.mutable_read_alleles())[StrCat(sample_id, "_alt_read_", i)] =
        alt_allele;

    Allele* new_allele =
        (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
    *new_allele = alt_allele;
  }

  return allele_count;
}

AlleleCount MakeTestMultiAlleleCount(
                int total_n,  // total number of reads
                const std::string& sample_id,  // sample id
                const std::string& ref,  // reference allele
                const absl::flat_hash_map<
                    std::string,
                    std::vector<std::string>>& alts_and_reads,  // map of alt
                                                // allele to supporting reads
                int start  // start position
                ) {
  AlleleCount allele_count;
  *(allele_count.mutable_position()) = nucleus::MakePosition("chr1", start);
  allele_count.set_ref_base(ref.substr(0, 1));
  int total_alt_n = 0;
  for (const auto& [alt, reads] : alts_and_reads) {
    total_alt_n += reads.size();
    AlleleType allele_type = AlleleTypeFromAlt(ref, alt);
    std::string alt_bases = alt;
    if (allele_type == AlleleType::DELETION) {
      alt_bases = ref.substr(0, ref.size() - alt.size() + 1);
    }
    if (allele_type == AlleleType::SUBSTITUTION ||
        allele_type == AlleleType::REFERENCE) {
      alt_bases = alt.substr(0, 1);
    }
    const Allele read_allele = MakeAllele(alt_bases,
                                          AlleleTypeFromAlt(ref, alt), 1);
    for (const auto& read : reads) {
      (*allele_count.mutable_read_alleles())[StrCat(sample_id, read)] =
          read_allele;

      Allele* new_allele =
          (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
      *new_allele = read_allele;
    }
  }
  CHECK_GE(total_n, total_alt_n) <<
      "Total number of reads must be >= n alt reads";
  allele_count.set_ref_supporting_read_count(total_n - total_alt_n);
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

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  VariantCallerOptions options = BasicOptions();
  int window_size = 5;
  int skip_next_count = 0;
  int prev_deletion_end = 0;
  options.set_small_model_vaf_context_window_size(window_size);
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options,
          target_sample_allele_counts,
          "sample");

  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] =
      caller->allele_counters_per_sample_["sample"]->Counts().begin() + 3;
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                         &allele_counter_iterators["sample"],
                         skip_next_count, prev_deletion_end);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
  EXPECT_THAT(call.allele_frequency_at_position().at(8), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(9), 85);
  EXPECT_THAT(call.allele_frequency_at_position().at(10), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(11), 100);
  EXPECT_THAT(call.allele_frequency_at_position().at(12), 0);
  caller->Clear();
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

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  int window_size = 3;
  VariantCallerOptions options = BasicOptions();
  options.set_small_model_vaf_context_window_size(window_size);

  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options,
          target_sample_allele_counts,
          "sample");
  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] =
      caller->allele_counters_per_sample_["sample"]->Counts().begin() + 3;

  int skip_next_count = 0;
  int prev_deletion_end = 0;
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                         &allele_counter_iterators["sample"],
                        skip_next_count, prev_deletion_end);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
  EXPECT_THAT(call.allele_frequency_at_position().at(9), 85);
  EXPECT_THAT(call.allele_frequency_at_position().at(10), 95);
  EXPECT_THAT(call.allele_frequency_at_position().at(11), 100);
  caller->Clear();
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

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {};
  allele_counts["sample"] = allele_count;

  int window_size = 0;
  VariantCallerOptions options = BasicOptions();
  options.set_small_model_vaf_context_window_size(window_size);

  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options,
          target_sample_allele_counts,
          "sample");
  // start in the middle of the iterator (3rd position);
  allele_counter_iterators["sample"] =
      caller->allele_counters_per_sample_["sample"]->Counts().begin() + 3;
  int skip_next_count = 0;
  int prev_deletion_end = 0;
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                         &allele_counter_iterators["sample"],
                        skip_next_count, prev_deletion_end);
  EXPECT_TRUE(static_cast<bool>(optional_variant));

  DeepVariantCall call = optional_variant.value();
  EXPECT_THAT(call.allele_frequency_at_position().size(), window_size);
  caller->Clear();
}

// Test AddSupportingReads with multiple samples.
TEST(VariantCallingTest,
     TestAddSupportingReadsWithMultipleSamples) {
  AlleleCount allele_count_1 =
      MakeTestAlleleCount(20, 19, "sample_1", "A", "T", 0);
  AlleleCount allele_count_2 =
      MakeTestAlleleCount(20, 5, "sample_2", "A", "C", 0);
  const std::vector<AlleleCount>& target_sample_allele_counts = {
      allele_count_1,
      allele_count_2,
  };
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;

  absl::node_hash_map<std::string, AlleleCount> allele_counts = {
    {"sample_1", allele_count_1},
    {"sample_2", allele_count_2},
  };

  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          BasicOptions(),
          target_sample_allele_counts,
          "sample_1");
  int skip_next_count = 0;
  int prev_deletion_end = 0;
  const std::optional<DeepVariantCall> optional_variant_1 =
      caller->CallVariant(allele_counts,
                         &allele_counter_iterators["sample_1"],
                         skip_next_count, prev_deletion_end);
  EXPECT_TRUE(static_cast<bool>(optional_variant_1));

  DeepVariantCall call = optional_variant_1.value();
  DeepVariantCall::SupportingReadsExt support_ext_1 =
      call.allele_support_ext().at("T");
  EXPECT_THAT(call.allele_support_ext().contains("C"), false);
  EXPECT_THAT(support_ext_1.read_infos().size(), 19);
  EXPECT_THAT(support_ext_1.read_infos().at(0).sample_name(), "sample_1");
  EXPECT_THAT(support_ext_1.read_infos().at(0).mapping_quality(), 30);
  EXPECT_THAT(support_ext_1.read_infos().at(0).average_base_quality(), 30);
  caller->Clear();
}

// Compute methylation statistics when only reference alleles are present.
TEST(VariantCallingTest, TestCallVariantComputeMethylationStatsOnlyReference) {
  VariantCallerOptions options = BasicOptions();
  // Currently, methylation calculations need ref
  // reads be tracked otherwise only alt reads are considered in
  // allele_counts.read_alleles.
  bool track_ref_reads = true;

  // Create an AlleleCount with 20 total reads, 0 alternate,
  // 2 methylated reference reads, and 0 methylated alternate read.
  AlleleCount allele_count = MakeTestAlleleCount(
      20, 10, "sample", "G", "C", 500, 2, 0, track_ref_reads);

  // Create target sample allele counts for context
  const std::vector<AlleleCount> target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "G", "C", 497),
      MakeTestAlleleCount(20, 7, "sample", "G", "C", 498),
      MakeTestAlleleCount(20, 9, "sample", "G", "C", 499),
      allele_count,  // Position of interest
      MakeTestAlleleCount(20, 11, "sample", "G", "C", 501),
      MakeTestAlleleCount(20, 13, "sample", "G", "C", 502),
      MakeTestAlleleCount(20, 15, "sample", "G", "C", 503),
  };

  // Create VariantCaller
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options, target_sample_allele_counts, "sample");

  absl::node_hash_map<std::string, AlleleCount> allele_counts;
  allele_counts["sample"] = allele_count;

  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;

  int skip_next_count = 0;
  int prev_deletion_end = 0;

  // Compute Methylation Stats via CallVariant
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                          &allele_counter_iterators["sample"],
                          skip_next_count, prev_deletion_end);

  DeepVariantCall call = optional_variant.value();

  // // Ensure variant has expected reference and alternate bases
  EXPECT_EQ(call.variant().reference_bases(), "G");

  // Ensure INFO field contains MF and MD
  EXPECT_TRUE(call.variant().calls(0).info().contains("MF"));
  EXPECT_TRUE(call.variant().calls(0).info().contains("MD"));

  // Check MF (Methylation Fraction)
  // MF = (# methylated ref) / (total ref = total - alt)
  // MF (REF) = (2) / (20-10) = 0.2
  EXPECT_EQ(call.variant().calls(0).info().at("MF").values(0).number_value(),
            0.2);
  // MF (ALT) = (0) / (10) = 0.0
  EXPECT_EQ(call.variant().calls(0).info().at("MF").values(1).number_value(),
            0);

  // Check MD (Methylation Depth)
  // MD = number of methylated reads
  // MD = 2 (methylated ref)
  EXPECT_EQ(call.variant().calls(0).info().at("MD").values(0).int_value(), 2);
  // MD = 0 (methylated alt)
  EXPECT_EQ(call.variant().calls(0).info().at("MD").values(1).int_value(), 0);

  caller->Clear();
}

// Calculate the methylation statistics when some alleles are methylated.
TEST(VariantCallingTest, TestCallVariantComputeMethylationStatsSomeMethylated) {
  VariantCallerOptions options = BasicOptions();
  // Currently, methylation calculations need ref
  // reads be tracked otherwise only alt reads are considered in
  // allele_counts.read_alleles.
  bool track_ref_reads = true;

  // Create an AlleleCount with 20 total reads, 10 alternate,
  // 2 methylated reference reads, and 3 methylated alternate reads.
  AlleleCount allele_count = MakeTestAlleleCount(
      20, 10, "sample", "G", "C", 500, 2, 3, track_ref_reads);

  // Create target sample allele counts for context
  const std::vector<AlleleCount> target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "G", "C", 497),
      MakeTestAlleleCount(20, 7, "sample", "G", "C", 498),
      MakeTestAlleleCount(20, 9, "sample", "G", "C", 499),
      allele_count,  // Position of interest
      MakeTestAlleleCount(20, 11, "sample", "G", "C", 501),
      MakeTestAlleleCount(20, 13, "sample", "G", "C", 502),
      MakeTestAlleleCount(20, 15, "sample", "G", "C", 503),
  };

  // Create VariantCaller
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options, target_sample_allele_counts, "sample");

  absl::node_hash_map<std::string, AlleleCount> allele_counts;
  allele_counts["sample"] = allele_count;

  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;

  int skip_next_count = 0;
  int prev_deletion_end = 0;

  // Compute Methylation Stats via CallVariant
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                          &allele_counter_iterators["sample"],
                          skip_next_count, prev_deletion_end);

  DeepVariantCall call = optional_variant.value();

  // Ensure variant has expected reference and alternate bases
  EXPECT_EQ(call.variant().reference_bases(), "G");

  // Ensure INFO field contains MF and MD
  EXPECT_TRUE(call.variant().calls(0).info().contains("MF"));
  EXPECT_TRUE(call.variant().calls(0).info().contains("MD"));

  // Check MF (Methylation Fraction)
  // MF (REF) = (# methylated ref) / (total ref)
  // MF (REF) = (2) / (10) = 0.2
  EXPECT_EQ(call.variant().calls(0).info().at("MF").values(0).number_value(),
            0.2);

  // MF (ALT) = (# methylated alt) / (total alt)
  // MF (ALT) = (3) / (10) = 0.3
  EXPECT_EQ(call.variant().calls(0).info().at("MF").values(1).number_value(),
            0.3);

  // Check MD (Methylation Depth)
  // MD (REF) = 2 methylated ref reads
  EXPECT_EQ(call.variant().calls(0).info().at("MD").values(0).int_value(), 2);

  // MD (ALT) = 3 methylated alt reads
  EXPECT_EQ(call.variant().calls(0).info().at("MD").values(1).int_value(), 3);

  caller->Clear();
}

// Calculate the methylation statistics when all alleles are not methylated.
TEST(VariantCallingTest, TestCallVariantComputeMethylationStatsNoMethylation) {
  VariantCallerOptions options = BasicOptions();
  // Currently, methylation calculations need ref
  // reads be tracked otherwise only alt reads are considered in
  // allele_counts.read_alleles.
  bool track_ref_reads = true;

  // Create an AlleleCount with 20 total reads, 10 alternate,
  // 0 methylated reference reads, and 0 methylated alternate reads.
  AlleleCount allele_count = MakeTestAlleleCount(
      20, 10, "sample", "G", "C", 500, 0, 0, track_ref_reads);

  // Create target sample allele counts for context
  const std::vector<AlleleCount> target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "G", "C", 497),
      MakeTestAlleleCount(20, 7, "sample", "G", "C", 498),
      MakeTestAlleleCount(20, 9, "sample", "G", "C", 499),
      allele_count,  // Position of interest
      MakeTestAlleleCount(20, 11, "sample", "G", "C", 501),
      MakeTestAlleleCount(20, 13, "sample", "G", "C", 502),
      MakeTestAlleleCount(20, 15, "sample", "G", "C", 503),
  };

  // Create VariantCaller
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options, target_sample_allele_counts, "sample");

  absl::node_hash_map<std::string, AlleleCount> allele_counts;
  allele_counts["sample"] = allele_count;

  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;

  int skip_next_count = 0;
  int prev_deletion_end = 0;

  // Compute Methylation Stats via CallVariant
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                          &allele_counter_iterators["sample"],
                          skip_next_count, prev_deletion_end);

  DeepVariantCall call = optional_variant.value();

  // Ensure variant has expected reference and alternate bases
  EXPECT_EQ(call.variant().reference_bases(), "G");

  // Ensure INFO field contains MF and MD
  EXPECT_FALSE(call.variant().calls(0).info().contains("MF"));
  EXPECT_FALSE(call.variant().calls(0).info().contains("MD"));

  caller->Clear();
}

// Ensure low-quality methylated ALT reads are excluded from methylation counts.
TEST(VariantCallingTest,
     TestCallVariantComputeMethylationStatsIgnoreLowQualityMethylation) {
  VariantCallerOptions options = BasicOptions();
  bool track_ref_reads = true;

  // Create an AlleleCount with 20 total reads, 10 alternate,
  // 3 methylated alternate reads, and 2 of them are low-quality.
  AlleleCount allele_count = MakeTestAlleleCount(
      20, 10, "sample", "G", "C", 500, 0, 3, track_ref_reads);
  allele_count.mutable_read_alleles()
      ->at("sample_alt_read_1")
      .set_is_low_quality(true);
  allele_count.mutable_read_alleles()
      ->at("sample_alt_read_2")
      .set_is_low_quality(true);

  // Create target sample allele counts for context
  const std::vector<AlleleCount> target_sample_allele_counts = {
      MakeTestAlleleCount(20, 5, "sample", "G", "C", 497),
      MakeTestAlleleCount(20, 7, "sample", "G", "C", 498),
      MakeTestAlleleCount(20, 9, "sample", "G", "C", 499),
      allele_count,  // Position of interest
      MakeTestAlleleCount(20, 11, "sample", "G", "C", 501),
      MakeTestAlleleCount(20, 13, "sample", "G", "C", 502),
      MakeTestAlleleCount(20, 15, "sample", "G", "C", 503),
  };

  // Create VariantCaller
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options, target_sample_allele_counts, "sample");

  absl::node_hash_map<std::string, AlleleCount> allele_counts;
  allele_counts["sample"] = allele_count;

  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;

  int skip_next_count = 0;
  int prev_deletion_end = 0;

  // Compute Methylation Stats via CallVariant
  const std::optional<DeepVariantCall> optional_variant =
      caller->CallVariant(allele_counts,
                          &allele_counter_iterators["sample"],
                          skip_next_count, prev_deletion_end);

  DeepVariantCall call = optional_variant.value();

  // Ensure variant has expected reference and alternate bases
  EXPECT_EQ(call.variant().reference_bases(), "G");

  // Ensure INFO field contains MF and MD with correct values
  EXPECT_TRUE(call.variant().calls(0).info().contains("MF"));
  EXPECT_TRUE(call.variant().calls(0).info().contains("MD"));

  nucleus::genomics::v1::ListValue expected_list_value;

  // Extract MF values and check correctness
  const auto& mf_list = call.variant().calls(0).info().at("MF").values();
  ASSERT_EQ(mf_list.size(), 2);
  EXPECT_FLOAT_EQ(mf_list[0].number_value(), 0.0);
  EXPECT_FLOAT_EQ(mf_list[1].number_value(), 0.125);

  // Extract MD values and check correctness
  const auto& md_list = call.variant().calls(0).info().at("MD").values();
  ASSERT_EQ(md_list.size(), 2);
  EXPECT_EQ(md_list[0].int_value(), 0);
  EXPECT_EQ(md_list[1].int_value(), 1);

  caller->Clear();
}

struct CreateComplexAllelesSupportTestData {
  absl::flat_hash_map<std::string, std::vector<AlleleAtPosition>>
          read_to_alt_alleles;
  int del_start;
  int del_len;
  std::string del_allele_ref_bases;
  absl::flat_hash_map<std::string, std::vector<std::string>>
      expected_complex_alleles_support;
};

class CreateComplexAllelesSupportTest
    : public testing::TestWithParam<CreateComplexAllelesSupportTestData> {};

TEST_P(CreateComplexAllelesSupportTest, CreateComplexAllelesSupportTestCases) {
  VariantCallerOptions options = BasicOptions();
  const VariantCaller caller(options);
  const CreateComplexAllelesSupportTestData& param = GetParam();
  auto complex_alleles_support = CreateComplexAllelesSupport(
      param.read_to_alt_alleles, param.del_start, param.del_len,
      param.del_allele_ref_bases);
  for (const auto& [complex_allele, reads] : complex_alleles_support) {
    auto it = param.expected_complex_alleles_support.find(complex_allele);
    ASSERT_NE(it, param.expected_complex_alleles_support.end());
    EXPECT_THAT(reads, testing::UnorderedElementsAreArray(it->second));
  }
  EXPECT_EQ(complex_alleles_support.size(),
            param.expected_complex_alleles_support.size());
}

INSTANTIATE_TEST_SUITE_P(
    CreateComplexAllelesSupportTests, CreateComplexAllelesSupportTest,
    testing::ValuesIn(std::vector<CreateComplexAllelesSupportTestData>({
      {
        .read_to_alt_alleles = {
          {"read1", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read2", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read3", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"T", AlleleType::REFERENCE, 13}}}
        },
        .del_start = 8,
        .del_len = 6,
        .del_allele_ref_bases = "CCGAATG",
        .expected_complex_alleles_support = {{"CCAAACG", {"read1", "read2"}},
                                             {"CCAAATG", {"read3"}}}
      },
      {
        .read_to_alt_alleles = {
          {"read1", {{"GATT", AlleleType::INSERTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read2", {{"GATT", AlleleType::INSERTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read3", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"T", AlleleType::REFERENCE, 13}}}
        },
        .del_start = 8,
        .del_len = 6,
        .del_allele_ref_bases = "CCGAATG",
        .expected_complex_alleles_support = {{"CCGATTAACG", {"read1", "read2"}},
                                             {"CCAAATG", {"read3"}}}
      },
      {
         // allele each.
        .read_to_alt_alleles = {
          {"read1", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read2", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"A", AlleleType::SUBSTITUTION, 13}}},
          {"read3", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"T", AlleleType::REFERENCE, 13}}}
        },
        .del_start = 8,
        .del_len = 6,
        .del_allele_ref_bases = "CCGAATG",
        .expected_complex_alleles_support = {{"CCAAACG", {"read1"}},
                                             {"CCAAAAG", {"read2"}},
                                             {"CCAAATG", {"read3"}}}
      },
      {
        .read_to_alt_alleles = {
          {"read1", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"C", AlleleType::SUBSTITUTION, 15}}},
          {"read2", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"A", AlleleType::SUBSTITUTION, 13}}},
          {"read3", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"T", AlleleType::REFERENCE, 13}}}
        },
        .del_start = 8,
        .del_len = 6,
        .del_allele_ref_bases = "CCGAATG",
        .expected_complex_alleles_support = {}
      },
      {
        .read_to_alt_alleles = {
          {"read1", {{"A", AlleleType::SUBSTITUTION, 8},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read2", {{"A", AlleleType::SUBSTITUTION, 8},
                     {"C", AlleleType::SUBSTITUTION, 13}}},
          {"read3", {{"A", AlleleType::SUBSTITUTION, 10},
                     {"T", AlleleType::REFERENCE, 13}}}
        },
        .del_start = 8,
        .del_len = 6,
        .del_allele_ref_bases = "CCGAATG",
        .expected_complex_alleles_support = {{"ACGAACG", {"read1", "read2"}},
                                             {"CCAAATG", {"read3"}}}
      }
    })));


struct ComplexVariantTestData {
  std::vector<AlleleCount> allele_count_context;
  int allele_counts_index = 0;
  int prev_deletion_end = 0;
  std::vector<Allele> expected_alleles;
};

class ComplexVariantTest
    : public testing::TestWithParam<ComplexVariantTestData> {};

TEST_P(ComplexVariantTest, ComplexVariantTestCases) {
  VariantCallerOptions options = BasicOptions();
  const ComplexVariantTestData& param = GetParam();
  std::unique_ptr<VariantCaller> caller =
      VariantCaller::MakeTestVariantCallerFromAlleleCounts(
          options, param.allele_count_context, "sample");
  std::string ref_bases;
  SelectAltAllelesResult ret = caller->SelectAltAlleles(
      {
        .allele_counts_by_sample = {
        {"sample",
         param.allele_count_context[param.allele_counts_index]}},
        .create_complex_alleles = true,
        .prev_deletion_end = param.prev_deletion_end,
      });
  EXPECT_THAT(ret.alt_alleles,
              testing::UnorderedPointwise(nucleus::EqualsProto(),
                                                   param.expected_alleles));

  // TODO: Also check modified allele_counts from retutrn value.
  caller->Clear();
}

INSTANTIATE_TEST_SUITE_P(
    ComplexVariantTests, ComplexVariantTest,
    testing::ValuesIn(std::vector<ComplexVariantTestData>({
      // Test case 1: Deletion allele overlaps one SNP.
      {
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTGGATCA",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "ACTGGATCA",
                  {"read_4", "read_5", "read_6"}
                }
              }, 7),
        MakeTestMultiAlleleCount(
          20, "sample", "ACTGGATCA",  // Ref allele
          {
          {}
          }, 8),
        MakeTestMultiAlleleCount(
          20, "sample", "T",  // SNP
          {
            {
              "G",
              {"read_4", "read_5", "read_6"}
            }
          }, 9)
      },
      .expected_alleles = {
          MakeAllele("ACTGGATCA", AlleleType::DELETION, 3),
          MakeAllele("ACGGGATCA", AlleleType::SUBSTITUTION, 3)
      }
      },
    { // Test case 2: Deletion allele overlaps two SNP.
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTGGATCA",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "ACTGGATCA",
                  {"read_4", "read_5", "read_6"}
                }
              }, 7
          ),
        MakeTestMultiAlleleCount(
          20, "sample", "T",  // SNP
          {
            {
              "G",
              {"read_4", "read_5", "read_6"}
            }
          }, 9),
        MakeTestMultiAlleleCount(
          20, "sample", "A",  // SNP
          {
            {
              "T",
              {"read_4", "read_5", "read_6"}
            }
          }, 12)
      },
      .expected_alleles = {
          MakeAllele("ACTGGATCA", AlleleType::DELETION, 3),
          MakeAllele("ACGGGTTCA", AlleleType::SUBSTITUTION, 3)
      }
    }
    // Test case 3: Deletion allele overlaps two SNPs supported by different
    // set of reads.
    // Ref:NNNNNNNAC T GG A TCANNNNNNN
    // read_1: NNNA- - -- - ---NNNNNNN
    // read_2: NNNA- - -- - ---NNNNNNN
    // read_3: NNNA- - -- - ---NNNNNNN
    // read_4: NNNAC[G]GG A TCANNNNNNN
    // read_5: NNNAC[G]GG A TCANNNNNNN
    // read_6: NNNAC[G]GG A TCANNNNNNN
    // read_7: NNNAC T GG[T]TCANNNNNNN
    // read_8: NNNAC T GG[T]TCANNNNNNN
    // read_9: NNNAC T GG[T]TCANNNNNNN
    , {
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTGGATCA",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "ACTGGATCA",
                  {"read_4", "read_5", "read_6", "read_7", "read_8", "read_9"}
                }
              }, 7)
        , MakeTestMultiAlleleCount(
          20, "sample", "T",  // SNP
          {
            {
              "G",
              {"read_4", "read_5", "read_6"}
            }
          }, 9)
        , MakeTestMultiAlleleCount(
          20, "sample", "A",  // SNP
          {
            {
              "T",
              {"read_7", "read_8", "read_9"}
            }
          }, 12)
      },
      .expected_alleles = {
          MakeAllele("ACTGGATCA", AlleleType::DELETION, 3),
          MakeAllele("ACGGGATCA", AlleleType::SUBSTITUTION, 3),
          MakeAllele("ACTGGTTCA", AlleleType::SUBSTITUTION, 3)
      }
    }
    , {  // Test case 4: Deletion allele overlaps SNP and insertion.
      // Ref:NNNNNNNAC T GG A --TCANNNNNNN
      // read_1: NNNA- - -- -   ---NNNNNNN
      // read_2: NNNA- - -- -   ---NNNNNNN
      // read_3: NNNA- - -- -   ---NNNNNNN
      // read_4: NNNAC[G]GG A   TCANNNNNNN
      // read_5: NNNAC[G]GG A   TCANNNNNNN
      // read_6: NNNAC[G]GG A   TCANNNNNNN
      // read_7: NNNAC T GG[ATT]TCANNNNNNN
      // read_8: NNNAC T GG[ATT]TCANNNNNNN
      // read_9: NNNAC T GG[ATT]TCANNNNNNN
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTGGATCA",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "ACTGGATCA",
                  {"read_4", "read_5", "read_6", "read_7", "read_8", "read_9"}
                }
              }, 7)
        , MakeTestMultiAlleleCount(
          20, "sample", "T",  // SNP
          {
            {
              "G",
              {"read_4", "read_5", "read_6"}
            }
          }, 9)
        , MakeTestMultiAlleleCount(
          20, "sample", "A",  // INS
          {
            {
              "ATT",
              {"read_7", "read_8", "read_9"}
            }
          }, 12)
      },
      .expected_alleles = {
          MakeAllele("ACTGGATCA", AlleleType::DELETION, 3),
          MakeAllele("ACGGGATCA", AlleleType::SUBSTITUTION, 3),
          MakeAllele("ACTGGATTTCA", AlleleType::SUBSTITUTION, 3)
      }
    }
    , {  // Test case 5: One base deletion allele and SNP right before the
    // deletion. This is not a complex variant, it is handled as a normal
    // multi-allele candidate.
      // Ref:NNNNNNN A CTGGATCANNNNNNN
      // read_1: NNN A C-GGATCANNNNNNN
      // read_2: NNN A C-GGATCANNNNNNN
      // read_3: NNN A C-GGATCANNNNNNN
      // read_4: NNN[T]CTGGATCANNNNNNN
      // read_5: NNN[T]CTGGATCANNNNNNN
      // read_6: NNN[T]CTGGATCANNNNNNN
      .allele_count_context = {
        MakeTestMultiAlleleCount(
          20, "sample", "A",  // SNP
          {
            {
              "T",
              {"read_4", "read_5", "read_6"}
            }
          }, 7),
        MakeTestMultiAlleleCount(
              20, "sample", "CT",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
              }, 8)
      },
      .allele_counts_index = 1,
      .expected_alleles = {
          MakeAllele("CT", AlleleType::DELETION, 3),
      }
    }
    , {  // Test case 6: One base deletion allele overlaps with SNP.
      // Ref:NNNNNNNA C TGGATCANNNNNNN
      // read_1: NNNA - TGGATCANNNNNNN
      // read_2: NNNA - TGGATCANNNNNNN
      // read_3: NNNA - TGGATCANNNNNNN
      // read_4: NNNA[A]TGGATCANNNNNNN
      // read_5: NNNA[A]TGGATCANNNNNNN
      // read_6: NNNA[A]TGGATCANNNNNNN
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "AC",  // Del allele
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "AC",
                  {"read_4", "read_5", "read_6"}
                }
              }, 7)
        , MakeTestMultiAlleleCount(
          20, "sample", "C",  // SNP
          {
            {
              "A",
              {"read_4", "read_5", "read_6"}
            }
          }, 8)
      },
      .expected_alleles = {
          MakeAllele("AC", AlleleType::DELETION, 3),
          MakeAllele("AA", AlleleType::SUBSTITUTION, 3),
      }
    }
    , {  // Test case 7: Multiple deletions of different lengths starting at
    // the same position.
      // Ref:NNNNNNNACTGGATCANNNNNNN
      // read_1: NNNA---GATCANNNNNNN
      // read_2: NNNA---GATCANNNNNNN
      // read_3: NNNA---GATCANNNNNNN
      // read_4: NNNA--GGATCANNNNNNN
      // read_5: NNNA--GGATCANNNNNNN
      // read_6: NNNA--GGATCANNNNNNN
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTG",  // 3 base deletion
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "AG",
                  {"read_4", "read_5", "read_6"}
                },
              }, 7)
      },
      .expected_alleles = {
          MakeAllele("ACTG", AlleleType::DELETION, 3),
          MakeAllele("ACT", AlleleType::DELETION, 3),
      }
    }
    , {  // Test case 8: Multiple deletions of different lengths starting at
    // the same position + SNP. In that case complex variant is not created.
    // Alleles created for position 7 will only contain two deletions and no
    // SNP.
      // Ref:NNNNNNNACTGGATCANNNNNNN
      // read_1: NNNA---GATCANNNNNNN
      // read_2: NNNA---GATCANNNNNNN
      // read_3: NNNA---GATCANNNNNNN
      // read_4: NNNA--GGATCANNNNNNN
      // read_5: NNNA--GGATCANNNNNNN
      // read_6: NNNA--GGATCANNNNNNN
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTG",  // 3 base deletion
              {
                {
                  "A",
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "AG",
                  {"read_4", "read_5", "read_6"}
                },
              }, 7)
        , MakeTestMultiAlleleCount(
            20, "sample", "C", {}, 8  // Ref
        )
        , MakeTestMultiAlleleCount(
            20, "sample", "T", {}, 9  // Ref
        )
        , MakeTestMultiAlleleCount(
          20, "sample", "G",  // SNP
          {
            {
              "A",
              {"read_4", "read_5", "read_6"}
            }
          }, 10)
      },
      .expected_alleles = {
          MakeAllele("ACTG", AlleleType::DELETION, 3),
          MakeAllele("ACT", AlleleType::DELETION, 3),
      }
    },
    {  // Test case 9: Complex variant is not created if previous deletion
    // overlaps current position. This test case has the same input as test
    // case 1 with the exception that prev_deletion_end is set to 8.
      .allele_count_context = {
        MakeTestMultiAlleleCount(
              20, "sample", "ACTGGATCA",
              {
                {
                  "A",  // Del allele
                  {"read_1", "read_2", "read_3"}
                },
                {
                  "ACTGGATCA",  // Ref allele
                  {"read_4", "read_5", "read_6"}
                }
              }, 7),
        MakeTestMultiAlleleCount(
          20, "sample", "ACTGGATCA",  // Ref allele
          {
          {}
          }, 8),
        MakeTestMultiAlleleCount(
          20, "sample", "T",  // SNP
          {
            {
              "G",
              {"read_4", "read_5", "read_6"}
            }
          }, 9)
      },
      .prev_deletion_end = 8,
      .expected_alleles = {
          MakeAllele("ACTGGATCA", AlleleType::DELETION, 3)
      }
    },
    })));

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
