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

#include "deepvariant/variant_calling.h"

#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/strings/string_view.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "google/protobuf/repeated_field.h"
#include "google/protobuf/text_format.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace vcf_candidate_importer {

using nucleus::EqualsProto;
using nucleus::MakePosition;
using nucleus::MakeRange;
using nucleus::VcfReader;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using ::testing::DoubleNear;
using ::testing::Eq;
using ::testing::UnorderedElementsAre;

constexpr char kSampleName[] = "MySampleName";
constexpr char kChr[] = "chr1";
constexpr int64_t kStart = 10;

AlleleCount MakeAlleleCount(const absl::string_view chr_name, int start,
                            const absl::string_view ref_base,
                            int ref_supporting_read_count,
                            const std::vector<Allele>& read_alleles) {
  AlleleCount allele_count;
  *allele_count.mutable_position() = nucleus::MakePosition(chr_name, start);
  QCHECK_EQ(ref_base.length(), 1) << "AlleleCount.ref_base has to be one char.";
  allele_count.set_ref_base(ref_base.data(), ref_base.size());
  allele_count.set_ref_supporting_read_count(ref_supporting_read_count);
  for (int i = 0; i < read_alleles.size(); ++i) {
    (*allele_count.mutable_read_alleles())[absl::StrCat("read_", i)] =
        read_alleles[i];
  }
  return allele_count;
}

AlleleCount MakeTestAlleleCount(int total_n, int alt_n,
                                const absl::string_view ref = "A",
                                int start = 100) {
  CHECK_GE(total_n, alt_n) << "Total number of reads must be >= n alt reads";
  std::vector<Allele> read_alleles;
  const Allele read_allele = MakeAllele("C", AlleleType::SUBSTITUTION, 1);
  for (int i = 0; i < alt_n; ++i) {
    read_alleles.push_back(read_allele);
  }
  return MakeAlleleCount(kChr, start, ref, total_n - alt_n, read_alleles);
}

enum class ExpectedVariant {
  kNoVariantExpected,
  kVariantExpected,
  kMaybeExpected,
};

VariantCallerOptions MakeOptions(
    const int min_count = 0, const double min_fraction = 0.0,
    absl::string_view sample_name = kSampleName,
    const double fraction_reference_sites_to_emit = -1.0) {
  VariantCallerOptions options;
  options.set_min_count_snps(min_count);
  options.set_min_count_indels(min_count);
  options.set_min_fraction_snps(min_fraction);
  options.set_min_fraction_indels(min_fraction);
  options.set_sample_name(std::string(sample_name));
  if (fraction_reference_sites_to_emit > 0)
    options.set_fraction_reference_sites_to_emit(
        fraction_reference_sites_to_emit);
  options.set_p_error(0.01);
  options.set_max_gq(50);
  options.set_gq_resolution(1);
  options.set_ploidy(2);

  return options;
}

Variant MakeExpectedVariant(const absl::string_view ref,
                            const std::vector<absl::string_view>& alts,
                            const int64_t start = kStart) {
  Variant variant;
  variant.set_reference_name(kChr);
  variant.set_start(start);
  variant.set_reference_bases(ref.data(), ref.size());
  for (const auto alt_allele : alts)
    variant.add_alternate_bases(alt_allele.data(), alt_allele.size());

  if (alts.empty()) {
    // Variant should be a single bp gVCF record with the kGVCFAltAllele
    // marker, genotypes of 0/0, and a GQ value of 0 (currently not
    // determined). Note that these are simple, baseline tests for any site
    // that doesn't have a variant call. Detailed testing of the proper gVCF
    // statistical calculations will come when those calculations appear.
    // TODO: Revise this test in full gVCF calculations CL.
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
  return MakeExpectedVariant(ref, {});
}

class VariantCallingTest : public ::testing::Test {
 protected:
  void CheckCall(const absl::string_view ref, const int min_alt_count,
                 const std::vector<Allele>& alleles,
                 const ExpectedVariant expect_variant,
                 const Variant& partial_expected_variant) {
    CheckCall(ref, VariantCaller(MakeOptions(min_alt_count)), alleles,
              expect_variant, partial_expected_variant);
  }

  // TODO: The CheckCallFromComputeVariantWithSameProposedVariant
  // functions below should be cleaned up. And, the assumption that the
  // proposed variant is the same as the expected variant should also be
  // re-examined.
  std::optional<DeepVariantCall>
  CheckCallFromComputeVariantWithSameProposedVariant(
      absl::string_view ref, const VariantCaller& caller,
      const std::vector<Allele>& alleles, const ExpectedVariant expect_variant,
      const Variant& expected_variant) {
    AlleleCount allele_count = ConstructAlleleCount(ref, alleles);
    std::vector<AlleleCount> allele_counts = {allele_count};
    const Variant& proposed_variant = expected_variant;
    return CheckCallFromComputeVariant(caller, proposed_variant, allele_counts,
                                       expect_variant, expected_variant);
  }

  void CheckCallFromComputeVariantWithSameProposedVariant(
      absl::string_view ref, const int min_alt_count,
      const std::vector<Allele>& alleles, const ExpectedVariant expect_variant,
      const Variant& partial_expected_variant) {
    CheckCallFromComputeVariantWithSameProposedVariant(
        ref, VariantCaller(MakeOptions(min_alt_count)), alleles, expect_variant,
        partial_expected_variant);
  }

  // Checks the result of CallVariant on an AlleleCount with the requested
  // properties from the arguments. Returns the resulting DeepVariantCall
  // produced by CallVariants for further testing in the callee.
  std::optional<DeepVariantCall> CheckCall(const absl::string_view ref,
                                      const VariantCaller& caller,
                                      const std::vector<Allele>& alleles,
                                      const ExpectedVariant expect_variant,
                                      const Variant& expected_variant) {
    AlleleCount allele_count = ConstructAlleleCount(ref, alleles);
    const std::optional<DeepVariantCall> optional_variant =
        caller.CallVariant(allele_count);
    CheckVariant(optional_variant, expect_variant, expected_variant);
    return optional_variant;
  }

  // Checks the result of ComputeVariant on an AlleleCount, Variant with the
  // requested properties from the arguments. Returns the resulting
  // DeepVariantCall produced by ComputeVariants for further testing in the
  // callee.
  std::optional<DeepVariantCall> CheckCallFromComputeVariant(
      const VariantCaller& caller, const Variant& proposed_variant,
      const std::vector<AlleleCount>& allele_counts,
      const ExpectedVariant expect_variant, const Variant& expected_variant) {
    const std::optional<DeepVariantCall> optional_variant =
        caller.ComputeVariant(proposed_variant, allele_counts);
    CheckVariant(optional_variant, expect_variant, expected_variant);
    return optional_variant;
  }

  AlleleCount ConstructAlleleCount(absl::string_view ref,
                                   const std::vector<Allele>& alleles) {
    // Construct the synthetic AlleleCount we'll use to call.
    AlleleCount allele_count;
    *allele_count.mutable_position() = MakePosition(kChr, kStart);
    allele_count.set_ref_base(std::string(ref));
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
          const string read_name = absl::StrCat("read_", ++read_counter);
          (*allele_count.mutable_read_alleles())[read_name] = read_allele;
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
    for (const absl::string_view ref : {"A", "C", "G", "T"}) {
      CheckCall(ref, 3, {MakeAllele(ref, AlleleType::REFERENCE, count)},
                ExpectedVariant::kNoVariantExpected, NoVariant(ref));
    }
  }
}

TEST_F(VariantCallingTest, TestNoVariantFromSoftclips) {
  for (const int count : {0, 1, 10, 100}) {
    const Allele allele = MakeAllele("ACCCCC", AlleleType::SOFT_CLIP, count);
    CheckCall("A", 3, {allele}, ExpectedVariant::kNoVariantExpected,
              NoVariant());
  }
}

TEST_F(VariantCallingTest, TestSNP) {
  for (const int count : {10, 100}) {
    for (const absl::string_view ref : {"A", "C", "G", "T"}) {
      for (const absl::string_view alt : {"A", "C", "G", "T"}) {
        if (alt != ref) {
          const Variant variant = MakeExpectedVariant(ref, {alt});
          // there's just alt observed
          CheckCall(ref, 3, {MakeAllele(alt, AlleleType::SUBSTITUTION, count)},
                    ExpectedVariant::kVariantExpected,
                    WithCounts(variant, {0, count}));
          // we see ref and alt, result is still the same
          CheckCall(ref, 3,
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
  const string alt = "C";
  const Allele alt_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, count);

  // If ref is "A", a canonical base, we make a call.
  CheckCall("A", 3, {alt_allele}, ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant("A", {alt}), {0, count}));
  // If ref isn't a canonical base, we don't make a call.
  CheckCall("N", 3, {alt_allele}, ExpectedVariant::kNoVariantExpected,
            NoVariant());
  CheckCall("R", 3, {alt_allele}, ExpectedVariant::kNoVariantExpected,
            NoVariant());
}

TEST_F(VariantCallingTest, TestMinCount1) {
  const int count = 10;
  const string ref = "A";
  const string alt = "C";
  const Variant variant =
      WithCounts(MakeExpectedVariant(ref, {alt}), {0, count});

  const Allele alt_allele = MakeAllele(alt, AlleleType::SUBSTITUTION, count);
  CheckCall(ref, count + 1, {alt_allele}, ExpectedVariant::kNoVariantExpected,
            NoVariant());
  CheckCall(ref, count, {alt_allele}, ExpectedVariant::kVariantExpected,
            variant);
  CheckCall(ref, count - 1, {alt_allele}, ExpectedVariant::kVariantExpected,
            variant);
}

TEST_F(VariantCallingTest, TestMinCount2) {
  const int count = 10;
  const string ref = "A";
  const string alt1 = "C";
  const string alt2 = "G";

  // Alt1 is above threshold and alt2 is below, so our variant has only alt1.
  CheckCall(
      ref, count,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count - 1),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1}), {0, count}, 2 * count - 1));

  // Both alt1 and alt2 are above threshold, so we get a multi-allelic back.
  CheckCall(
      ref, count,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}), {0, count, count}));

  // Both alt1 and alt2 are below the threshold, so we get a no-call
  CheckCall(ref, count,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count - 1),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count - 1),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
}

TEST_F(VariantCallingTest, TestMinFraction1) {
  const int count = 10;
  const string ref = "A";
  const string alt = "C";
  const Variant variant = MakeExpectedVariant(ref, {alt});

  VariantCaller caller(MakeOptions(count, 0.1));

  // just an alt, with count >= the min so it goes in
  CheckCall(ref, caller, {MakeAllele(alt, AlleleType::SUBSTITUTION, count)},
            ExpectedVariant::kVariantExpected, WithCounts(variant, {0, count}));

  // both ref and alt are above the min count and are at 0.5 fraction >= 0.1
  // so our record is variant
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {count, count}));

  // Here ref is so frequent that alt goes below our 0.1 fraction so the
  // result is a no-call
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count * 100),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());

  // Checking the symmetric case : alt is very frequent
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count),
                MakeAllele(alt, AlleleType::SUBSTITUTION, count * 100),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {count, 100 * count}));
}

TEST_F(VariantCallingTest, TestMinFractionMultiAllelic) {
  const int count = 10;
  const string ref = "A";
  const string alt1 = "C";
  const string alt2 = "G";
  VariantCaller caller(MakeOptions(count, 0.1));

  // both alt1 and alt2 are above the min count and are at 0.5 fraction >= 0.1
  // so our record is contains both
  CheckCall(
      ref, caller,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}), {0, count, count}));

  // Here alt1 is so frequent that alt2 goes below our 0.1 fraction so the
  // result is a no-call
  CheckCall(ref, caller,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count * 100),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant(ref, {alt1}), {0, count * 100},
                       count * 101));

  // Checking the symmetric case : alt2 is very frequent
  CheckCall(ref, caller,
            {
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count * 100),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(MakeExpectedVariant(ref, {alt2}), {0, count * 100},
                       count * 101));

  // Finally, ref is so frequent that neither alt1 or alt2 are good enough
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, count * 100),
                MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
                MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
}

TEST_F(VariantCallingTest, TestMinSNPIndelSeparately) {
  const string ref = "A";
  const string snp_alt = "C";
  const string ins_alt = "AC";
  const string del_alt = "AC";
  const Variant snp_variant = MakeExpectedVariant(ref, {snp_alt});
  const Variant ins_variant = MakeExpectedVariant(ref, {ins_alt});
  const Variant del_variant = MakeExpectedVariant(del_alt, {ref});

  VariantCallerOptions options;
  options.set_min_count_snps(5);
  options.set_min_count_indels(10);
  options.set_min_fraction_snps(0.1);
  options.set_min_fraction_indels(0.5);
  options.set_sample_name(kSampleName);
  options.set_ploidy(2);
  VariantCaller caller(options);

  // Check that we respect min_count for SNPs and indels. With a count of 8 we
  // satisfy the requirements for SNPs but not for indels. Bumping the indel
  // allele counts to their minimum produces calls.
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(snp_alt, AlleleType::SUBSTITUTION, 8),
            },
            ExpectedVariant::kVariantExpected, WithCounts(snp_variant, {8, 8}));
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(ins_alt, AlleleType::INSERTION, 8),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(ins_alt, AlleleType::INSERTION,
                           options.min_count_indels()),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(ins_variant, {8, options.min_count_indels()}));
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 8),
                MakeAllele(del_alt, AlleleType::DELETION, 8),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(
      ref, caller,
      {
          MakeAllele(ref, AlleleType::REFERENCE, 8),
          MakeAllele(del_alt, AlleleType::DELETION, options.min_count_indels()),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(del_variant, {8, options.min_count_indels()}));

  // Check that we respect min_fraction for SNPs and indels. With 20% of the
  // pileup having the allele we satisfy the requirements for SNPs but not for
  // indels. At 50% we call both the SNPs and the indels.
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(snp_alt, AlleleType::SUBSTITUTION, 20),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(snp_variant, {80, 20}));
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(ins_alt, AlleleType::INSERTION, 20),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(ins_alt, AlleleType::INSERTION, 80),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(ins_variant, {80, 80}));
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(del_alt, AlleleType::DELETION, 20),
            },
            ExpectedVariant::kNoVariantExpected, NoVariant());
  CheckCall(ref, caller,
            {
                MakeAllele(ref, AlleleType::REFERENCE, 80),
                MakeAllele(del_alt, AlleleType::DELETION, 80),
            },
            ExpectedVariant::kVariantExpected,
            WithCounts(del_variant, {80, 80}));
}

TEST_F(VariantCallingTest, TestMultAllelicSNP) {
  const int count = 10;
  const string ref = "A";
  const string alt1 = "C";
  const string alt2 = "G";
  const Variant variant = MakeExpectedVariant(ref, {alt1, alt2});
  CheckCall(ref, count,
            {MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
             MakeAllele(alt2, AlleleType::SUBSTITUTION, count)},
            ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count, count}));
}

TEST_F(VariantCallingTest, TestBiAllelicDeletion) {
  for (const absl::string_view alt_bases : {"AC", "ACCC", "ACCCCCCCCC"}) {
    const int count = 10;
    const string ref = "A";
    const Allele alt = MakeAllele(alt_bases, AlleleType::DELETION, count);
    const Variant variant = MakeExpectedVariant(alt.bases(), {ref});
    CheckCall(ref, count, {alt}, ExpectedVariant::kVariantExpected,
              WithCounts(variant, {0, count}));
  }
}

TEST_F(VariantCallingTest, TestBiAllelicInsertion) {
  for (const absl::string_view alt_bases : {"AC", "ACCC", "ACCCCCCCCC"}) {
    const int count = 10;
    const string ref = "A";
    const Allele alt = MakeAllele(alt_bases, AlleleType::INSERTION, count);
    const Variant variant = MakeExpectedVariant(ref, {alt.bases()});
    CheckCall(ref, count, {alt}, ExpectedVariant::kVariantExpected,
              WithCounts(variant, {0, count}));
  }
}

TEST_F(VariantCallingTest, TestDeletionInsertion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("ACCC", AlleleType::INSERTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "ACCCTGC"});
  CheckCall("A", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestTwoDeletions) {
  const int count = 10;
  const Allele alt1 = MakeAllele("AT", AlleleType::DELETION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "AGC"});
  CheckCall("A", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestTwoInsertions) {
  const int count = 10;
  const Allele alt1 = MakeAllele("AT", AlleleType::INSERTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::INSERTION, count + 1);
  const Variant variant = MakeExpectedVariant("A", {"AT", "ATGC"});
  CheckCall("A", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count, count + 1}));
}

TEST_F(VariantCallingTest, TestSNPDeletion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("C", AlleleType::SUBSTITUTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("ATGC", {"A", "CTGC"});
  CheckCall("A", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count + 1, count}));
}

TEST_F(VariantCallingTest, TestDeletionWithNonRefAnchor) {
  // In this case we have a deletion with a non-reference anchor base, which
  // produces a complex variant. Check that the conversion works correctly.
  // See http://internal for more information about this issue.
  const int count = 10;
  const Allele alt = MakeAllele("AA", AlleleType::DELETION, count);
  const Variant variant = MakeExpectedVariant("TA", {"A"});
  CheckCall("T", count, {alt}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count}));
}

TEST_F(VariantCallingTest, TestInsertionWithNonRefAnchor) {
  const int count = 10;
  const Allele alt = MakeAllele("AA", AlleleType::INSERTION, count);
  const Variant variant = MakeExpectedVariant("T", {"AA"});
  CheckCall("T", count, {alt}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count}));
}

TEST_F(VariantCallingTest, TestDeletionWithNonRefAnchor2) {
  // In this case we have two deletions, one with a non-reference anchor base,
  // which produces a complex variant. Check that the conversion works
  // correctly. See http://internal for more information about this issue.
  const int count = 10;
  const Allele alt1 = MakeAllele("AA", AlleleType::DELETION, count);
  const Allele alt2 = MakeAllele("TA", AlleleType::DELETION, count + 1);
  const Variant variant = MakeExpectedVariant("TA", {"A", "T"});
  CheckCall("T", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
            WithCounts(variant, {0, count, count + 1}));
}

TEST_F(VariantCallingTest, TestSNPInsertion) {
  const int count = 10;
  const Allele alt1 = MakeAllele("C", AlleleType::SUBSTITUTION, count);
  const Allele alt2 = MakeAllele("ATGC", AlleleType::INSERTION, count + 1);
  const Variant variant = MakeExpectedVariant("A", {"ATGC", "C"});
  CheckCall("A", count, {alt1, alt2}, ExpectedVariant::kVariantExpected,
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
  CheckCall("A", count, {alt1, alt2, alt3, alt4, alt5},
            ExpectedVariant::kVariantExpected,
            WithCounts(variant,
                       {0, count + 3, count + 1, count + 2, count + 4, count}));
}

// Extracts the read_names from the map value of call.allele_support at key,
// returning them as a vector of strings.
std::vector<std::string> SupportingReadNames(const DeepVariantCall& call,
                                             absl::string_view key) {
  std::vector<std::string> names;
  for (const string& read_name :
       call.allele_support().at(std::string(key)).read_names()) {
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
  const string ref = "A";
  const string alt1 = "ACT";
  const string alt2 = "ATG";
  const string alt3 = "G";
  const Variant variant = MakeExpectedVariant("ATG", {"A", "ACTTG"});
  VariantCaller caller(MakeOptions(count, 0.1));

  const std::optional<DeepVariantCall> optional_call =
      CheckCall(ref, caller,
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
  const string ref = "A";
  const string alt = "C";
  const Variant variant = MakeExpectedVariant("A", {"."});
  VariantCaller caller(MakeOptions(count, 0.1, kSampleName, 1.0));

  const std::optional<DeepVariantCall> optional_call = CheckCall(
      ref, caller,
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
  const string ref = "A";
  const Variant variant = MakeExpectedVariant("A", {"."});
  VariantCaller caller(MakeOptions(count, 0.1, kSampleName, fraction));

  const int tries = 10000;
  int successes = 0;
  for (int i = 0; i < tries; ++i) {
    const std::optional<DeepVariantCall> optional_call =
        CheckCall(ref, caller, {MakeAllele(ref, AlleleType::REFERENCE, count)},
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
  std::vector<AlleleCount> allele_counts = {
      MakeTestAlleleCount(0, 0, "A", 10), MakeTestAlleleCount(10, 10, "G", 11),
      MakeTestAlleleCount(0, 0, "G", 12), MakeTestAlleleCount(0, 0, "G", 13),
      MakeTestAlleleCount(11, 9, "T", 14)};

  const VariantCaller caller(MakeOptions());
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromAlleleCounts(allele_counts);

  // We expect our candidates to have 2 and 5 in order.
  Variant variant2 = WithCounts(MakeExpectedVariant("G", {"C"}, 11), {0, 10});
  Variant variant5 = WithCounts(MakeExpectedVariant("T", {"C"}, 14), {2, 9});
  ASSERT_THAT(candidates.size(), Eq(2));
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant2));
  EXPECT_THAT(candidates[1].variant(), EqualsProto(variant5));
}

TEST_F(VariantCallingTest, TestCallPositionsFromVcfQueryingVcf) {
  const VariantCaller caller(MakeOptions());
  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(
          nucleus::GetTestData("vcf_candidate_importer.indels.chr20.vcf.gz",
                               "deepvariant/testdata/input"),
          nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  std::vector<AlleleCount> allele_count_not_used = {AlleleCount()};

  // Querying in contigInHeaderWithCandidates returns one candidate.
  std::vector<int> positions = caller.CallPositionsFromVcf(
      allele_count_not_used, MakeRange("chr20", 59777020, 59974170),
      reader.get());
  EXPECT_THAT(positions, UnorderedElementsAre(
    59777552,
    59804672,
    59848583,
    59858359,
    59858388,
    59865297,
    59884411,
    59904401,
    59904404,
    59906637,
    59912353,
    59928658,
    59951038,
    59958677,
    59958833,
    59965315,
    59965632,
    59965720,
    59974165
    ));
  }

TEST_F(VariantCallingTest, TestCallsFromVcfQueryingVcf) {
  // TODO
  // For human readability, here is the content of the VCF:
  /*
  $ zcat learning/genomics/deepvariant/testdata/input/test_calls_from_vcf.vcf.gz
  ##fileformat=VCFv4.2
  ##contig=<ID=contigInHeaderWithCandidates,length=10>
  ##contig=<ID=contigInHeaderNoCandidates,length=10>
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
  contigInHeaderWithCandidates    3       .       T       G       60      PASS    .       GT      ./.
  contigNotInHeaderWithCandidates 1       .       A       G       60      PASS    .       GT      ./.
  */
  const VariantCaller caller(MakeOptions());
  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(
          nucleus::GetTestData("test_calls_from_vcf.vcf.gz",
                               "deepvariant/testdata/input"),
          nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  std::vector<AlleleCount> allele_count_not_used = {AlleleCount()};

  // Querying in contigInHeaderWithCandidates returns one candidate.
  std::vector<DeepVariantCall> candidates1 = caller.CallsFromVcf(
      allele_count_not_used, MakeRange("contigInHeaderWithCandidates", 0, 5),
      reader.get());
  EXPECT_EQ(candidates1.size(), 1);

  // Querying contigInHeaderNoCandidates returns 0 candidate and doesn't crash.
  std::vector<DeepVariantCall> candidates2 = caller.CallsFromVcf(
      allele_count_not_used, MakeRange("contigInHeaderNoCandidates", 0, 5),
      reader.get());
  EXPECT_EQ(candidates2.size(), 0);

  // Querying contigNotInVcf returns empty results.
  std::vector<DeepVariantCall> candidates3 = caller.CallsFromVcf(
      allele_count_not_used, MakeRange("contigNotInVcf", 0, 5), reader.get());
  EXPECT_EQ(candidates3.size(), 0);
}

TEST_F(VariantCallingTest, TestCallsFromVcfDetails) {
  // TODO
  // For human readability, here is the content of the VCF:
  /*
  $ zcat learning/genomics/deepvariant/testdata/input/test_calls_from_vcf.vcf.gz
  ##fileformat=VCFv4.2
  ##contig=<ID=contigInHeaderWithCandidates,length=10>
  ##contig=<ID=contigInHeaderNoCandidates,length=10>
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
  contigInHeaderWithCandidates    3       .       T       G       60      PASS    .       GT      ./.
  contigNotInHeaderWithCandidates 1       .       A       G       60      PASS    .       GT      ./.
  */
  const VariantCaller caller(MakeOptions());
  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(
          nucleus::GetTestData("test_calls_from_vcf.vcf.gz",
                               "deepvariant/testdata/input"),
          nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  // The AlleleCount represents 10 reads.
  // 5 reads supporting ref "T"
  // 3 reads supporting alt "A"
  // 2 reads supporting alt "G"
  // ---------------
  //   POS: 0123456
  //   REF:   T
  // ---------------
  // Reads:   A
  //          A
  //          A
  //          G
  //          G
  //          T
  //          T
  //          T
  //          T
  //          T
  // ---------------
  const AlleleCount allele_count =
      MakeAlleleCount("contigInHeaderWithCandidates",  // chr_name
                      2,                               // start
                      "T",                             // ref_base
                      5,  // ref_supporting_read_count
                      {MakeAllele("A", AlleleType::SUBSTITUTION, 1),
                       MakeAllele("A", AlleleType::SUBSTITUTION, 1),
                       MakeAllele("A", AlleleType::SUBSTITUTION, 1),
                       MakeAllele("G", AlleleType::SUBSTITUTION, 1),
                       MakeAllele("G", AlleleType::SUBSTITUTION, 1)});

  // Querying contigInHeaderWithCandidates returns 1 candidate.
  std::vector<DeepVariantCall> candidates = caller.CallsFromVcf(
      {allele_count}, MakeRange("contigInHeaderWithCandidates", 0, 5),
      reader.get());
  EXPECT_EQ(candidates.size(), 1);
  EXPECT_EQ(candidates[0].variant().reference_bases(), "T");
  EXPECT_EQ(candidates[0].variant().alternate_bases_size(), 1);
  EXPECT_EQ(candidates[0].variant().alternate_bases(0), "G");
  EXPECT_EQ(candidates[0].variant().calls_size(), 1);
  EXPECT_EQ(candidates[0].variant().calls(0).info_size(), 3);
  nucleus::genomics::v1::ListValue expected_list_value;
  google::protobuf::TextFormat::ParseFromString(
      "values { int_value: 5 } values { int_value: 2 }", &expected_list_value);
  EXPECT_THAT(candidates[0].variant().calls(0).info().at("AD"),
              EqualsProto(expected_list_value));
  google::protobuf::TextFormat::ParseFromString("values { int_value: 10 }",
                                      &expected_list_value);
  EXPECT_THAT(candidates[0].variant().calls(0).info().at("DP"),
              EqualsProto(expected_list_value));
  google::protobuf::TextFormat::ParseFromString("values { number_value: 0.2 }",
                                      &expected_list_value);
  EXPECT_THAT(candidates[0].variant().calls(0).info().at("VAF"),
              EqualsProto(expected_list_value));
}

TEST_F(VariantCallingTest, TestTrainUncalledGenotypes) {
  // For human readability, here is the content of the VCF:
  /*
  $ zcat learning/genomics/deepvariant/testdata/input/test_calls_from_vcf.vcf.gz
  ##fileformat=VCFv4.2
  ##contig=<ID=contigInHeaderWithCandidates,length=10>
  ##contig=<ID=contigInHeaderNoCandidates,length=10>
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  default
  contigInHeaderWithCandidates    3       .       T       G       60      PASS    .       GT      ./.
  contigNotInHeaderWithCandidates 1       .       A       G       60      PASS    .       GT      ./.
  */
  VariantCallerOptions options = MakeOptions();
  options.set_skip_uncalled_genotypes(true);
  const VariantCaller caller(options);
  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(
          nucleus::GetTestData("test_calls_from_vcf.vcf.gz",
                               "deepvariant/testdata/input"),
          nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  std::vector<AlleleCount> allele_count_not_used = {AlleleCount()};

  // Querying in contigInHeaderWithCandidates returns 0 candidates.
  std::vector<DeepVariantCall> candidates1 = caller.CallsFromVcf(
      allele_count_not_used, MakeRange("contigInHeaderWithCandidates", 0, 5),
      reader.get());
  EXPECT_EQ(candidates1.size(), 0);
}

TEST_F(VariantCallingTest, TestCallsFromVariantsInRegion) {
  // Our test AlleleCounts are 5 positions:
  //
  // 1: A ref  [no reads]
  // 2: G/C variant
  // 3: G ref  [no reads]
  // 4: G ref  [no reads]
  // 5: T/C variant
  //
  std::vector<AlleleCount> allele_counts = {
      MakeTestAlleleCount(0, 0, "A", 10), MakeTestAlleleCount(10, 10, "G", 11),
      MakeTestAlleleCount(0, 0, "G", 12), MakeTestAlleleCount(0, 0, "G", 13),
      MakeTestAlleleCount(11, 9, "T", 14)};

  Variant variant2 = WithCounts(MakeExpectedVariant("G", {"C"}, 11), {0, 10});
  Variant variant5 = WithCounts(MakeExpectedVariant("T", {"C"}, 14), {2, 9});
  std::vector<Variant> variants = {variant2, variant5};

  const VariantCaller caller(MakeOptions());
  std::vector<DeepVariantCall> candidates =
      caller.CallsFromVariantsInRegion(allele_counts, variants);

  // We expect our candidates to have 2 and 5 in order.
  ASSERT_THAT(candidates.size(), Eq(2));
  EXPECT_THAT(candidates[0].variant(), EqualsProto(variant2));
  EXPECT_THAT(
      SupportingReadNames(candidates[0], "C"),
      UnorderedElementsAre("read_1", "read_0", "read_2", "read_7", "read_8",
                           "read_6", "read_3", "read_4", "read_5", "read_9"));
  EXPECT_THAT(candidates[1].variant(), EqualsProto(variant5));
  EXPECT_THAT(
      SupportingReadNames(candidates[1], "C"),
      UnorderedElementsAre("read_1", "read_0", "read_2", "read_7", "read_8",
                           "read_6", "read_3", "read_4", "read_5"));
}

TEST_F(VariantCallingTest, TestComputeVariantMultiAllelic) {
  const int count = 10;
  const string ref = "A";
  const string alt1 = "C";
  const string alt2 = "G";
  VariantCaller caller(MakeOptions(count, 0.1));

  // Both alt1 and alt2 are above the min count and are at 0.5 fraction >= 0.1
  // so our record is contains both.
  CheckCallFromComputeVariantWithSameProposedVariant(
      ref, caller,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}), {0, count, count}));

  // Here alt1 is so frequent that alt2 goes below our 0.1 fraction so we only
  // call alt1.
  CheckCallFromComputeVariantWithSameProposedVariant(
      ref, caller,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count * 100),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1}), {0, count * 100},
                 count * 101));

  // Checking the symmetric case : alt2 is very frequent.
  CheckCallFromComputeVariantWithSameProposedVariant(
      ref, caller,
      {
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count * 100),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt2}), {0, count * 100},
                 count * 101));

  // Finally, ref is very frequent but alt1 and alt2 are still included.
  CheckCallFromComputeVariantWithSameProposedVariant(
      ref, caller,
      {
          MakeAllele(ref, AlleleType::REFERENCE, count * 100),
          MakeAllele(alt1, AlleleType::SUBSTITUTION, count * 50),
          MakeAllele(alt2, AlleleType::SUBSTITUTION, count * 50),
      },
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant(ref, {alt1, alt2}),
                 {count * 100, count * 50, count * 50}, count * 200));
}

// See internal.
TEST_F(VariantCallingTest, TestComputeVariantDifferentRefs) {
  int count = 2;
  VariantCaller caller(MakeOptions(count));
  int ref_supporting_read_count = 9;

  const std::vector<Allele> read_alleles = {
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CA", AlleleType::DELETION, 1),
      MakeAllele("CAA", AlleleType::DELETION, 1),
      MakeAllele("CAA", AlleleType::DELETION, 1),
      MakeAllele("CAA", AlleleType::DELETION, 1),
  };
  const AlleleCount allele_count =
      MakeAlleleCount(kChr,      // chr_name
                      92457968,  // start
                      "C",       // ref_base
                      ref_supporting_read_count, read_alleles);

  Variant proposed_variant = MakeExpectedVariant("CA", {"C"}, 92457968);
  CheckCallFromComputeVariant(
      caller, proposed_variant, {allele_count},
      ExpectedVariant::kVariantExpected,
      WithCounts(MakeExpectedVariant("CAA", {"CA"}, 92457968),
                 {ref_supporting_read_count, 6},
                 ref_supporting_read_count + read_alleles.size()));
}

// See internal.
TEST_F(VariantCallingTest, TestComputeVariantDifferentRefs2) {
  int count = 2;
  VariantCaller caller(MakeOptions(count));
  int ref_supporting_read_count = 8;

  const std::vector<Allele> read_alleles = {
      MakeAllele("TACAC", AlleleType::DELETION, 1),
      MakeAllele("TACAC", AlleleType::DELETION, 1),
      MakeAllele("TACAC", AlleleType::DELETION, 1),
      MakeAllele("TACAC", AlleleType::DELETION, 1),
  };
  const AlleleCount allele_count =
      MakeAlleleCount(kChr,      // chr_name
                      66618315,  // start
                      "T",       // ref_base
                      ref_supporting_read_count, read_alleles);

  Variant proposed_variant =
      MakeExpectedVariant("TACACACACAC", {"TACACAC", "T"}, 66618315);
  std::optional<DeepVariantCall> dv_call = CheckCallFromComputeVariant(
      caller, proposed_variant, {allele_count},
      ExpectedVariant::kVariantExpected,
      // Now, the 4 "TACAC" DELELTIONs above are correctly counted under
      // TACACACACAC->TACACAC.
      // This is because now we correct recognize that
      // TACAC->T
      // is equivalent to
      // TACACACACAC->TACACAC.
      WithCounts(MakeExpectedVariant("TACACACACAC", {"TACACAC", "T"}, 66618315),
                 {ref_supporting_read_count, 4, 0},
                 ref_supporting_read_count + read_alleles.size()));
  QCHECK_EQ(dv_call->allele_support_size(), 1);
  // Confirm that the 4 alleles "MakeAllele("TACAC", AlleleType::DELETION, 1)
  // above are correctly added to the corrresponding variant
  // "TACACACACAC->TACACAC", which is the same as "TACAC->T" after the
  // right-trimming simplification.
  const auto it = dv_call->allele_support().find("TACACAC");
  QCHECK(it != dv_call->allele_support().end());
  QCHECK_EQ(it->second.read_names_size(), 4);
}

}  // namespace vcf_candidate_importer
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
