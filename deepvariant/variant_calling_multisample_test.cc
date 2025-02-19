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

  // const VariantCaller caller(options);
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

  // const VariantCaller caller(options);
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
