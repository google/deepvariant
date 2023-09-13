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

#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/allelecounter.h"
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

using nucleus::EqualsProto;
using nucleus::MakePosition;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using ::testing::DoubleNear;
using ::testing::Eq;
using ::testing::UnorderedElementsAre;

constexpr char kSampleName[] = "MySampleName";
constexpr char kChr[] = "chr1";
constexpr int64_t kStart = 10;

// This function needs to return a pointer because it needs to be wrapped
// into nucleus::ConstProtoPtr. variant_calling_multisample API is only used by
// Python, we don't call CallsFromAlleleCounts from C++.
AlleleCount MakeTestAlleleCount(int total_n, int alt_n,
                                const std::string& sample_id,
                                const std::string& ref = "A",
                                const std::string& alt = "C", int start = 100,
                                const std::string& read_prefix = "") {
  CHECK_GE(total_n, alt_n) << "Total number of reads must be >= n alt reads";
  std::string test_read_prefix = read_prefix;
  if (test_read_prefix.empty()) {
    test_read_prefix = sample_id;
  }
  AlleleCount allele_count;
  *(allele_count.mutable_position()) = nucleus::MakePosition("chr1", start);
  allele_count.set_ref_base(ref);
  allele_count.set_ref_supporting_read_count(total_n - alt_n);
  const Allele read_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, 1);
  for (int i = 0; i < alt_n; ++i) {
    (*allele_count.mutable_read_alleles())[absl::StrCat(
        test_read_prefix, "_read_", i)] = read_allele;

    Allele* new_allele =
        (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
    *new_allele = read_allele;
  }
  return allele_count;
}

void ReleaseAlleleCounterPointers(
    const std::unordered_map<std::string, AlleleCounter*> allele_counters) {
  for (auto allele_counter_entry : allele_counters) {
    delete allele_counter_entry.second;
  }
}

enum class ExpectedVariant {
  kNoVariantExpected,
  kVariantExpected,
  kMaybeExpected,
};

VariantCallerOptions MakeOptions(
    const int min_count = 0, const double min_fraction = 0.0,
    const std::string& sample_name = kSampleName,
    const double fraction_reference_sites_to_emit = -1.0) {
  VariantCallerOptions options;
  options.set_min_count_snps(min_count);
  options.set_min_count_indels(min_count);
  options.set_min_fraction_snps(min_fraction);
  options.set_min_fraction_indels(min_fraction);
  options.set_min_fraction_multiplier(1.0);
  options.set_sample_name(sample_name);
  if (fraction_reference_sites_to_emit > 0)
    options.set_fraction_reference_sites_to_emit(
        fraction_reference_sites_to_emit);
  options.set_p_error(0.01);
  options.set_max_gq(50);
  options.set_gq_resolution(1);
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

// Creates a non-variant site with the given reference base (defaults to "A").
Variant NoVariant(absl::string_view ref = "A") {
  return MakeExpectedVariant(std::string(ref), {});
}

class VariantCallingTest : public ::testing::Test {
 protected:
  void CheckCall(const std::string& ref, const int expected_len,
                 const int min_alt_count, const std::vector<Allele>& alleles,
                 const ExpectedVariant expect_variant,
                 const Variant& partial_expected_variant) {
    CheckCall(ref, expected_len, VariantCaller(MakeOptions(min_alt_count)),
              alleles, expect_variant, partial_expected_variant);
  }

  // Checks the result of CallVariant on an AlleleCount with the requested
  // properties from the arguments. Returns the resulting DeepVariantCall
  // produced by CallVariants for further testing in the callee.
  std::optional<DeepVariantCall> CheckCall(const std::string& ref,
                                      const int expected_len,
                                      const VariantCaller& caller,
                                      const std::vector<Allele>& alleles,
                                      const ExpectedVariant expect_variant,
                                      const Variant& expected_variant) {
    AlleleCount allele_count = ConstructAlleleCount(ref, alleles, "sample_id");
    absl::node_hash_map<std::string, AlleleCount> allele_count_per_sample;
    allele_count_per_sample["sample_id"] = allele_count;
    const std::optional<DeepVariantCall> optional_variant =
        caller.CallVariant(allele_count_per_sample, "sample_id");
    CheckVariant(optional_variant, expect_variant, expected_variant);
    return optional_variant;
  }

  AlleleCount ConstructAlleleCount(const std::string& ref,
                                   const std::vector<Allele>& alleles,
                                   const std::string& sample_id) {
    // Construct the synthetic AlleleCount we'll use to call.
    AlleleCount allele_count;
    *allele_count.mutable_position() = MakePosition(kChr, kStart);
    allele_count.set_ref_base(ref);
    int read_counter = 0;
    for (const Allele& allele : alleles) {
      if (allele.type() == AlleleType::REFERENCE) {
        // Reference alleles are stored as counts in ref_supporting_read_count.
        int prev = allele_count.ref_supporting_read_count();
        allele_count.set_ref_supporting_read_count(prev + allele.count());
        // Ensure that we keep read names consistent across ref/non-ref reads.
        read_counter += allele.count();
      } else {
        // Non-reference reads are stored in the read_alleles list.
        const Allele read_allele = MakeAllele(allele.bases(), allele.type(), 1);
        for (int i = 0; i < allele.count(); ++i) {
          const std::string read_name = absl::StrCat("read_", ++read_counter);
          (*allele_count.mutable_read_alleles())[read_name] = read_allele;

          Allele* new_allele =
              (*allele_count.mutable_sample_alleles())[sample_id].add_alleles();
          *new_allele = allele;
        }
      }
    }
    return allele_count;
  }

  void CheckVariant(const std::optional<DeepVariantCall> optional_variant,
                    const ExpectedVariant expect_variant,
                    const Variant& partial_expected_variant) {
    switch (expect_variant) {
      case ExpectedVariant::kNoVariantExpected: {
        EXPECT_FALSE(static_cast<bool>(optional_variant));
        break;
      }
      case ExpectedVariant::kVariantExpected: {
        EXPECT_TRUE(static_cast<bool>(optional_variant));
        if (optional_variant) {
          // Checking optional_variant deals with our case where we really want
          // to ASSERT_THAT but ASSERT cannot be used in a helper with a
          // non-void return.
          EXPECT_THAT(optional_variant->variant(),
                      EqualsProto(partial_expected_variant));
        }
        break;
      }
      case ExpectedVariant::kMaybeExpected:
        // We may or may not make a call at this site, so there's nothing to
        // check.
        break;
      default:  // We don't have an expectation for any other options.
        LOG(FATAL) << "ExpectedVariant state: "
                   << static_cast<int>(expect_variant);
    }
  }
};

TEST_F(VariantCallingTest, TestNoVariant) {
  for (const int count : {0, 1, 10, 100}) {
    for (const std::string& ref : {"A", "C", "G", "T"}) {
      CheckCall(ref, 1, 3, {MakeAllele(ref, AlleleType::REFERENCE, count)},
                ExpectedVariant::kNoVariantExpected, NoVariant(ref));
    }
  }
}

TEST_F(VariantCallingTest, TestNoVariantFromSoftclips) {
  for (const int count : {0, 1, 10, 100}) {
    const Allele allele = MakeAllele("ACCCCC", AlleleType::SOFT_CLIP, count);
    CheckCall("A", 1, 3, {allele}, ExpectedVariant::kNoVariantExpected,
              NoVariant());
  }
}

TEST_F(VariantCallingTest, TestSNP) {
  for (const int count : {10, 100}) {
    for (const std::string& ref : {"A", "C", "G", "T"}) {
      for (const std::string& alt : {"A", "C", "G", "T"}) {
        if (alt != ref) {
          const Variant variant = MakeExpectedVariant(ref, {alt});
          // there's just alt observed
          CheckCall(ref, 1, 3,
                    {MakeAllele(alt, AlleleType::SUBSTITUTION, count)},
                    ExpectedVariant::kVariantExpected,
                    WithCounts(variant, {0, count}));
          // we see ref and alt, result is still the same
          CheckCall(ref, 1, 3,
                    {MakeAllele(alt, AlleleType::SUBSTITUTION, count),
                     MakeAllele(ref, AlleleType::REFERENCE, count)},
                    ExpectedVariant::kVariantExpected,
                    WithCounts(variant, {count, count}));
        }
      }
    }
  }
}

TEST_F(VariantCallingTest, TestNonCanonicalBase) {
  const int count = 100;
  const std::string alt = "C";
  const Allele alt_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, count);

  // If ref is "A", a canonical base, we make a call.
  CheckCall("A", 1, 3, {alt_allele}, ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant("A", {alt}), {0, count}));
  // If ref isn't a canonical base, we don't make a call.
  CheckCall("N", 1, 3, {alt_allele}, ExpectedVariant::kNoVariantExpected,
            NoVariant());
  CheckCall("R", 1, 3, {alt_allele}, ExpectedVariant::kNoVariantExpected,
            NoVariant());
}

TEST_F(VariantCallingTest, TestMinCount1) {
  const int count = 10;
  const std::string ref = "A";
  const std::string alt = "C";
  const Variant variant =
      WithCounts(MakeExpectedVariant(ref, {alt}), {0, count});

  const Allele alt_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, count);
  CheckCall(ref, 1, count + 1, {alt_allele},
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, 1, count, {alt_allele}, ExpectedVariant::kVariantExpected,
            variant);
  CheckCall(ref, 1, count - 1, {alt_allele}, ExpectedVariant::kVariantExpected,
            variant);
}

TEST_F(VariantCallingTest, TestMinCount2) {
  const int count = 10;
  const std::string ref = "A";
  const std::string alt1 = "C";
  const std::string alt2 = "G";

  // Alt1 is above threshold and alt2 is below, so our variant has only alt1.
  CheckCall(
      ref, 1, count,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count - 1),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1}), {0, count}, 2 * count - 1));

  // Both alt1 and alt2 are above threshold, so we get a multi-allelic back.
  CheckCall(
      ref, 1, count,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}), {0, count, count}));

  // Both alt1 and alt2 are below the threshold, so we get a no-call
  CheckCall(ref, 1, count,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count - 1),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count - 1),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
}

TEST_F(VariantCallingTest, TestMinFraction1) {
  const int count = 10;
  const std::string ref = "A";
  const std::string alt = "C";
  const Variant variant = MakeExpectedVariant(ref, {alt});

  VariantCaller caller(MakeOptions(count, 0.1));

  // just an alt, with count >= the min so it goes in
  CheckCall(ref, 1, caller, {MakeAllele(alt, AlleleType::SUBSTITUTION, count)},
            ExpectedVariant::kVariantExpected, WithCounts(variant, {0, count}));

  // both ref and alt are above the min count and are at 0.5 fraction >= 0.1
  // so our record is variant
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {count, count}));

  // Here ref is so frequent that alt goes below our 0.1 fraction so the
  // result is a no-call
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count * 100),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());

  // Checking the symmetric case : alt is very frequent
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count * 100),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {count, 100 * count}));
}

TEST_F(VariantCallingTest, TestMinFractionMultiAllelic) {
  const int count = 10;
  const std::string ref = "A";
  const std::string alt1 = "C";
  const std::string alt2 = "G";
  VariantCaller caller(MakeOptions(count, 0.1));

  // both alt1 and alt2 are above the min count and are at 0.5 fraction >= 0.1
  // so our record is contains both
  CheckCall(
      ref, 1, caller,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}), {0, count, count}));

  // Here alt1 is so frequent that alt2 goes below our 0.1 fraction so the
  // result is a no-call
  CheckCall(ref, 1, caller,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count * 100),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant(ref, {alt1}), {0, count * 100},
                       count * 101));

  // Checking the symmetric case : alt2 is very frequent
  CheckCall(ref, 1, caller,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count * 100),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant(ref, {alt2}), {0, count * 100},
                       count * 101));

  // Finally, ref is so frequent that neither alt1 or alt2 are good enough
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count * 100),
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
}

TEST_F(VariantCallingTest, TestMinSNPIndelSeparately) {
  const std::string ref = "A";
  const std::string snp_alt = "C";
  const std::string ins_alt = "AC";
  const std::string del_alt = "AC";
  const Variant snp_variant = MakeExpectedVariant(ref, {snp_alt});
  const Variant ins_variant = MakeExpectedVariant(ref, {ins_alt});
  const Variant del_variant = MakeExpectedVariant(del_alt, {ref});

  VariantCallerOptions options;
  options.set_min_count_snps(5);
  options.set_min_count_indels(10);
  options.set_min_fraction_snps(0.1);
  options.set_min_fraction_indels(0.5);
  options.set_min_fraction_multiplier(1.0);
  options.set_sample_name(kSampleName);
  options.set_ploidy(2);
  VariantCaller caller(options);

  // Check that we respect min_count for SNPs and indels. With a count of 8 we
  // satisfy the requirements for SNPs but not for indels. Bumping the indel
  // allele counts to their minimum produces calls.
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(snp_alt, AlleleType::SUBSTITUTION, 8),
            },
            ExpectedVariant::kVariantExpected, WithCounts(snp_variant, {8, 8}));
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(ins_alt, AlleleType::INSERTION, 8),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(ins_alt, AlleleType::INSERTION,
                           options.min_count_indels()),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(ins_variant, {8, options.min_count_indels()}));
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(del_alt, AlleleType::DELETION, 8),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(
      ref, 1, caller,
      {
          MakeAllele(ref, AlleleType::REFERENCE, 8),
          MakeAllele(del_alt, AlleleType::DELETION, options.min_count_indels()),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(del_variant, {8, options.min_count_indels()}));

  // Check that we respect min_fraction for SNPs and indels. With 20% of the
  // pileup having the allele we satisfy the requirements for SNPs but not for
  // indels. At 50% we call both the SNPs and the indels.
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(snp_alt, AlleleType::SUBSTITUTION, 20),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(snp_variant, {80, 20}));
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(ins_alt, AlleleType::INSERTION, 20),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(ins_alt, AlleleType::INSERTION, 80),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(ins_variant, {80, 80}));
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(del_alt, AlleleType::DELETION, 20),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, 1, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(del_alt, AlleleType::DELETION, 80),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(del_variant, {80, 80}));
}

TEST_F(VariantCallingTest, TestMultAllelicSNP) {
  const int count = 10;
  const std::string ref = "A";
  const std::string alt1 = "C";
  const std::string alt2 = "G";
  const Variant variant = MakeExpectedVariant(ref, {alt1, alt2});
  CheckCall(ref, 1, count,
            {MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
             MakeAllele(alt2, AlleleType::SUBSTITUTION, count)},
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count, count}));
}

TEST_F(VariantCallingTest, TestBiAllelicDeletion) {
  for (const std::string& alt_bases : {"AC", "ACCC", "ACCCCCCCCC"}) {
    const int count = 10;
    const std::string ref = "A";
    const Allele alt = MakeAllele(alt_bases, AlleleType::DELETION, count);
    const Variant variant = MakeExpectedVariant(alt.bases(), {ref});
    CheckCall(ref, alt_bases.size(), count, {alt},
              ExpectedVariant::kVariantExpected,
              WithCounts(variant, {0, count}));
  }
}

TEST_F(VariantCallingTest, TestBiAllelicInsertion) {
  for (const std::string& alt_bases : {"AC", "ACCC", "ACCCCCCCCC"}) {
    const int count = 10;
    const std::string ref = "A";
    const Allele alt = MakeAllele(alt_bases, AlleleType::INSERTION, count);
    const Variant variant = MakeExpectedVariant(ref, {alt.bases()});
    CheckCall(ref, 1, count, {alt}, ExpectedVariant::kVariantExpected,
              WithCounts(variant, {0, count}));
  }
}

TEST_F(VariantCallingTest, TestDeletionInsertion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("ACCC", AlleleType::INSERTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "ACCCTGC"});
  CheckCall("A", 4, count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestTwoDeletions) {
  const int count = 10;
  const Allele alt1 = MakeAllele("AT", AlleleType::DELETION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "AGC"});
  CheckCall("A", 4, count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestTwoInsertions) {
  const int count = 10;
  const Allele alt1 = MakeAllele("AT", AlleleType::INSERTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::INSERTION, count + 1);
  const Variant variant = MakeExpectedVariant("A", {"AT", "ATGC"});
  CheckCall("A", 1, count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count, count + 1}));
}

TEST_F(VariantCallingTest, TestSNPDeletion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("C", AlleleType::SUBSTITUTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "CTGC"});
  CheckCall("A", 4, count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestDeletionWithNonRefAnchor) {
  // In this case we have a deletion with a non-reference anchor base, which
  // produces a complex variant. Check that the conversion works correctly.
  // See http://internal for more information about this issue.
  const int count = 10;
  const Allele alt = MakeAllele("AA", AlleleType::DELETION, count);
  const Variant variant = MakeExpectedVariant("TA", {"A"});
  CheckCall("T", 4, count, {alt}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count}));
}

TEST_F(VariantCallingTest, TestInsertionWithNonRefAnchor) {
  const int count = 10;
  const Allele alt = MakeAllele("AA", AlleleType::INSERTION, count);
  const Variant variant = MakeExpectedVariant("T", {"AA"});
  CheckCall("T", 4, count, {alt}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count}));
}

TEST_F(VariantCallingTest, TestDeletionWithNonRefAnchor2) {
  // In this case we have two deletions, one with a non-reference anchor base,
  // which produces a complex variant. Check that the conversion works
  // correctly. See http://internal for more information about this issue.
  const int count = 10;
  const Allele alt1 = MakeAllele("AA", AlleleType::DELETION, count);
  const Allele alt2 = MakeAllele("TA", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("TA", {"T"});
  CheckCall("T", 4, count + 1, {alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1}));
}

// Following tests are not exported to the open source because I couldn't find
// an OSS replacement for util::tuple::testing::FieldPairsAre.

TEST_F(VariantCallingTest, TestSNPInsertion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("C", AlleleType::SUBSTITUTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::INSERTION, count + 1);
  const Variant variant = MakeExpectedVariant("A", {"ATGC", "C"});
  CheckCall("A", 1, count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestKitchenSink) {
  const int count = 10;
  const Allele alt1 = MakeAllele("C", AlleleType::SUBSTITUTION, count);
  const Allele alt2 = MakeAllele("AA", AlleleType::INSERTION, count + 1);
  const Allele alt3 = MakeAllele("ACAC", AlleleType::INSERTION, count + 2);
  const Allele alt4 = MakeAllele("ATGC", AlleleType::DELETION, count + 3);
  const Allele alt5 = MakeAllele("AT", AlleleType::DELETION, count + 4);
  const Variant variant = MakeExpectedVariant(
      "ATGC", {
                  // order changed due to sorting of alt alleles
                  "A",
                  "AATGC",
                  "ACACTGC",
                  "AGC",
                  "CTGC",
              });
  CheckCall("A", 4, count, {alt1, alt2, alt3, alt4, alt5},
            ExpectedVariant::kVariantExpected,
            WithCounts(variant,
                       {0, count + 3, count + 1, count + 2, count + 4, count}));
}

// Extracts the read_names from the map value of call.allele_support at key,
// returning them as a vector of strings.
std::vector<std::string> SupportingReadNames(const DeepVariantCall& call,
                                             const std::string& key) {
  std::vector<std::string> names;
  for (const std::string& read_name :
       call.allele_support().at(key).read_names()) {
    names.push_back(read_name);
  }
  return names;
}

TEST_F(VariantCallingTest, TestReadSupport) {
  // We have reads that support ref, alt1, alt2 and alt3 in the pileup. Alt3
  // doesn't have enough support to be a real alt allele. Because there are
  // insertion and deletion alleles we have a complex mapping between input
  // alleles and the resulting Variant alleles. Some reads support ref and won't
  // show up in the support map; the reads supporting alt3 get mapped to
  // supporting the kSupportingUncalledAllele allele, and the reads for the
  // insertion and deletion need to map properly from their initial read
  // alleles to different variant alleles.
  const int count = 5;
  const std::string ref = "A";
  const std::string alt1 = "ACT";
  const std::string alt2 = "ATG";
  const std::string alt3 = "G";
  const Variant variant = MakeExpectedVariant("ATG", {"A", "ACTTG"});
  VariantCaller caller(MakeOptions(count, 0.1));

  const std::optional<DeepVariantCall> optional_call =
      CheckCall(ref, 1, caller,
                {
                    MakeAllele(ref, AlleleType::REFERENCE, count),
                    MakeAllele(alt1, AlleleType::INSERTION, count),
                    MakeAllele(alt2, AlleleType::DELETION, count + 1),
                    MakeAllele(alt3, AlleleType::SUBSTITUTION, count - 1),
                },
                ExpectedVariant::kVariantExpected,
                WithCounts(variant, {count, count + 1, count}, 4 * count));
  ASSERT_TRUE(optional_call);
  const DeepVariantCall& call = *optional_call;

  // These inline read names implicitly know how the read names are generated in
  // CheckCall. Slightly ugly but allows us to very explicitly test the values
  // in the DeepVariantCall.allele_support map without doing any clever
  // calculations that are hard to do given the complex mapping between input
  // read alleles and output variant alleles.
  std::vector<std::string> keys;
  for (auto& entry : call.allele_support()) {
    keys.push_back(entry.first);
  }

  EXPECT_THAT(keys,
              UnorderedElementsAre("A", "ACTTG", kSupportingUncalledAllele));
  EXPECT_THAT(SupportingReadNames(call, "A"),
              UnorderedElementsAre("read_11", "read_12", "read_13", "read_14",
                                   "read_15", "read_16"));
  EXPECT_THAT(
      SupportingReadNames(call, "ACTTG"),
      UnorderedElementsAre("read_6", "read_7", "read_8", "read_9", "read_10"));
  EXPECT_THAT(SupportingReadNames(call, kSupportingUncalledAllele),
              UnorderedElementsAre("read_17", "read_18", "read_19", "read_20"));
}

TEST_F(VariantCallingTest, TestRefSites) {
  // Test that the reference variant call is coming back well formatted.
  const int count = 5;
  const std::string ref = "A";
  const std::string alt = "C";
  const Variant variant = MakeExpectedVariant("A", {"."});
  VariantCaller caller(MakeOptions(count, 0.1, kSampleName, 1.0));

  const std::optional<DeepVariantCall> optional_call = CheckCall(
      ref, 1, caller,
      {
          MakeAllele(ref, AlleleType::REFERENCE, count),
          MakeAllele(alt, AlleleType::SUBSTITUTION, 1),
      },
      ExpectedVariant::kVariantExpected, WithCounts(variant, {}, count + 1));
  ASSERT_TRUE(optional_call);

  const DeepVariantCall& call = *optional_call;
  EXPECT_THAT(SupportingReadNames(call, kSupportingUncalledAllele),
              UnorderedElementsAre(absl::StrCat("read_", count + 1)));
}

TEST_F(VariantCallingTest, TestRefSitesFraction) {
  // Test that the caller respects our requested fraction of ref sites.
  const double fraction = 0.6;
  const int count = 5;
  const std::string ref = "A";
  const Variant variant = MakeExpectedVariant("A", {"."});
  VariantCaller caller(MakeOptions(count, 0.1, kSampleName, fraction));

  const int tries = 10000;
  int successes = 0;
  for (int i = 0; i < tries; ++i) {
    const std::optional<DeepVariantCall> optional_call = CheckCall(
        ref, 1, caller, {MakeAllele(ref, AlleleType::REFERENCE, count)},
        ExpectedVariant::kMaybeExpected, variant);
    if (optional_call) {
      ++successes;
    }
  }

  const double rate = (1.0 * successes) / tries;
  EXPECT_THAT(rate, DoubleNear(fraction, 0.02));
}

TEST_F(VariantCallingTest, TestCallsFromAlleleCounts) {
  // Our test AlleleCounts are 5 positions:
  //
  // 1: A ref  [no reads]
  // 2: G/C variant
  // 3: G ref  [no reads]
  // 4: G ref  [no reads]
  // 5: T/C variant
  //
  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"sample_id",
       AlleleCounter::InitFromAlleleCounts(
           {MakeTestAlleleCount(0, 0, "sample_id", "A", "C", 10),
            MakeTestAlleleCount(10, 10, "sample_id", "G", "C", 11),
            MakeTestAlleleCount(0, 0, "sample_id", "G", "C", 12),
            MakeTestAlleleCount(0, 0, "sample_id", "G", "C", 13),
            MakeTestAlleleCount(11, 9, "sample_id", "T", "C", 14)})}};

  const VariantCaller caller(MakeOptions());
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "sample_id");

  // We expect our candidates to have 2 and 5 in order.
  Variant variant2 = WithCounts(MakeExpectedVariant("G", {"C"}, 11), {0, 10});
  Variant variant5 = WithCounts(MakeExpectedVariant("T", {"C"}, 14), {2, 9});
  ASSERT_THAT(candidates.size(), Eq(2));
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant2));
  EXPECT_THAT(candidates[1].variant(), EqualsProto(variant5));
  ReleaseAlleleCounterPointers(allele_counters);
}

// Testing that candidate is created for a target sample if ref support is very
// high in another sample. In which case allele fraction ratio would be too low
// for this candidate if we calculate allele ratio from all samples combined.
// This is to test the case when we call parent in DeepTrio and one of other
// samples have a very high coverage comparable to the target sample.
TEST_F(VariantCallingTest, TestCallsFromAlleleCountsUnevenCoverage) {
  // Parent 1:
  // SNP A -> T with 2 reads support and 2 reads ref support
  //
  // Child:
  // no alt alleles, 3 reads support ref
  //
  // Parent 2:
  // no alt alleles, 3 reads support ref
  //
  // set_min_fraction_snps > 2/(4 + 3 + 3)

  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"parent_1", AlleleCounter::InitFromAlleleCounts(
                       {MakeTestAlleleCount(4, 2, "parent_1", "A", "T", 10)})},
      {"child", AlleleCounter::InitFromAlleleCounts(
                    {MakeTestAlleleCount(3, 0, "child", "A", "A", 10)})},
      {"parent_2", AlleleCounter::InitFromAlleleCounts(
                       {MakeTestAlleleCount(3, 0, "parent_2", "A", "A", 10)})}};

  const VariantCaller caller(MakeOptions(0, 0.2));
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "parent_1");

  // We expect one candidate with 2 reads supporting alt and 2 reads supporting
  // ref.
  Variant variant = WithCounts(MakeExpectedVariant("A", {"T"}, 10), {2, 2});
  EXPECT_EQ(candidates.size(), 1);
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant));
  ReleaseAlleleCounterPointers(allele_counters);
}

// Testing that candidate is not created when allele fraction is less than
// threshold, and is created when allele_fraction_trio_coefficient is used.
TEST_F(VariantCallingTest,
       TestCallsFromAlleleCountsWithAlleleFractionTrioCoefficient) {
  // Parent 1:
  // SNP A -> T with 2 reads support and 18 reads ref support
  //
  // Child:
  // SNP A -> T with 1 reads support and 9 reads ref support
  //
  // Parent 2:
  // no alt alleles, 3 reads support ref
  //
  // set_min_fraction_snps >= 0.2

  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"parent_1", AlleleCounter::InitFromAlleleCounts(
                       {MakeTestAlleleCount(20, 2, "parent_1", "A", "T", 10)})},
      {"child", AlleleCounter::InitFromAlleleCounts(
                    {MakeTestAlleleCount(10, 7, "child", "A", "T", 10)})},
      {"parent_2", AlleleCounter::InitFromAlleleCounts({MakeTestAlleleCount(
                       10, 0, "parent_2", "A", "A", 10)})}};

  // When min_fraction_multiplier is set to 1.0 it doesn't affect
  // allele filter. We expect that candidate won't be created since allele
  // fraction is 0.225 is less than 0.226
  VariantCallerOptions options = MakeOptions(0, 0.226);
  options.set_min_fraction_multiplier(1.0);
  const VariantCaller caller(options);
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "parent_1");
  EXPECT_EQ(candidates.size(), 0);

  // Now we set min_fraction_multiplier to 0.5 which should result
  // in candidate to be created.
  VariantCallerOptions options2 = MakeOptions(0, 0.226);
  options2.set_min_fraction_multiplier(0.5);
  const VariantCaller caller2(options2);
  Variant variant = WithCounts(MakeExpectedVariant("A", {"T"}, 10), {18, 2});

  std::vector<DeepVariantCall> candidates2 =
      caller2.CallsFromAlleleCounts(allele_counters, "parent_1");
  EXPECT_EQ(candidates2.size(), 1);
  EXPECT_THAT(candidates2[0].variant(), EqualsProto(variant));

  ReleaseAlleleCounterPointers(allele_counters);
}

// Testing that candidate is created with low alt allele support if there is
// a support for the same allele in other samples.
TEST_F(VariantCallingTest,
       TestCallsFromAlleleCountsAltSupportInMultipleSamples) {
  // Parent 1:
  // SNP A -> T with 2 reads support and 18 reads ref support
  //
  // Child:
  // SNP A -> T with 1 reads support and 9 reads ref support
  //
  // Parent 2:
  // no alt alleles, 3 reads support ref
  //
  // set_min_fraction_snps >= 0.2

  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"parent_1", AlleleCounter::InitFromAlleleCounts(
                       {MakeTestAlleleCount(20, 2, "parent_1", "A", "T", 10)})},
      {
          "child",
          AlleleCounter::InitFromAlleleCounts(
              {MakeTestAlleleCount(10, 7, "child", "A", "T", 10)}),
      },
      {"parent_2", AlleleCounter::InitFromAlleleCounts({MakeTestAlleleCount(
                       10, 0, "parent_2", "A", "A", 10)})}};

  const VariantCaller caller(MakeOptions(0, 0.1999));
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "parent_1");

  // We expect our candidate to have 2 reads supporting alt and 18 reads
  // supporting ref.
  Variant variant = WithCounts(MakeExpectedVariant("A", {"T"}, 10), {18, 2});
  EXPECT_EQ(candidates.size(), 1);
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant));
  ReleaseAlleleCounterPointers(allele_counters);
}

// This test verifies that we can handle the case when one of the samples
// doesn't have any reads for some positions and as a result no AlleleCount
// objects are created for those positions. It is an overkill because
// empty AlleleCount is still created for those positions.
TEST_F(VariantCallingTest,
       TestCallsFromAlleleCountsAltSupportSomePositionsMissAlleleCounter) {
  // Position 10:
  // Parent 1:
  // SNP A -> T with 2 reads support alt and 18 reads ref support
  //
  // Child:
  // SNP A -> T with 1 reads support alt and 9 reads ref support
  //
  // Parent 2:
  // no alt alleles, 3 reads support ref
  //
  // set_min_fraction_snps >= 0.2
  //
  // Position 11:
  // Parent1:
  // no alt alleles, 3 reads alt support ref
  //
  // Child:
  // no alt alleles, 3 reads alt support ref
  //
  // Parent2:
  // no reads.

  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"parent_1", AlleleCounter::InitFromAlleleCounts({
                       MakeTestAlleleCount(20, 2, "parent_1", "A", "T", 10),
                       MakeTestAlleleCount(3, 0, "parent_1", "A", "A", 11),
                   })},
      {"child", AlleleCounter::InitFromAlleleCounts({
                    MakeTestAlleleCount(10, 7, "child", "A", "T", 10),
                    MakeTestAlleleCount(3, 0, "child", "A", "A", 11),
                })},
      {"parent_2", AlleleCounter::InitFromAlleleCounts({MakeTestAlleleCount(
                       10, 0, "parent_2", "A", "A", 10)})}};

  const VariantCaller caller(MakeOptions(0, 0.1999));
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "parent_1");

  // We expect our candidate to have 2 reads supporting alt and 18 reads
  // supporting ref.
  Variant variant = WithCounts(MakeExpectedVariant("A", {"T"}, 10), {18, 2});
  EXPECT_EQ(candidates.size(), 1);
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant));
  ReleaseAlleleCounterPointers(allele_counters);
}

// Testing that read with the same id coming from different samples is not
// used multiple times in read  support.
TEST_F(VariantCallingTest, TestCallsFromAlleleCountsDuplicateReadIds) {
  // Parent 1:
  // SNP A -> T with 3 reads support and 1 reads ref support
  //
  // Child:
  // SNP A -> T with 2 reads support and 1 reads ref support
  //
  // Parent 2:
  // no alt alleles, 3 reads support ref
  //

  const std::unordered_map<std::string, AlleleCounter*> allele_counters = {
      {"parent_1", AlleleCounter::InitFromAlleleCounts({MakeTestAlleleCount(
                       4, 3, "parent_1", "A", "T", 10, "common_prefix")})},
      {"child", AlleleCounter::InitFromAlleleCounts({MakeTestAlleleCount(
                    3, 2, "child", "A", "T", 10, "common_prefix")})},
      {"parent_2", AlleleCounter::InitFromAlleleCounts(
                       {MakeTestAlleleCount(3, 0, "parent_2", "A", "A", 10)})}};

  const VariantCaller caller(MakeOptions());
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counters, "child");

  // We expect one candidate with 3 reads (instead of 5) supporting alt.
  EXPECT_THAT(candidates, testing::SizeIs(1));
  EXPECT_THAT(
      SupportingReadNames(candidates[0], "T"),
      UnorderedElementsAre("common_prefix_read_0", "common_prefix_read_1",
                           "common_prefix_read_2"));
  ReleaseAlleleCounterPointers(allele_counters);
}

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
