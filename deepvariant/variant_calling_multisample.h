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

// A very simple but highly sensitive variant caller.
//
#ifndef LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_MULTISAMPLE_H_
#define LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_MULTISAMPLE_H_

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/check.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/samplers.h"

namespace nucleus {
class VcfReader;
}
namespace learning {
namespace genomics {
namespace deepvariant {
namespace multi_sample {

using learning::genomics::deepvariant::Allele;
using learning::genomics::deepvariant::AlleleCount;
using learning::genomics::deepvariant::AlleleCounter;
using learning::genomics::deepvariant::AlleleType;
using learning::genomics::deepvariant::DeepVariantCall;
using learning::genomics::deepvariant::VariantCallerOptions;
using nucleus::genomics::v1::Variant;

// The alternate allele string for the gVCF "any" alternate allele.
extern const char* const kGVCFAltAllele;

// In a DeepVariantCall, reads can support an allele that didn't pass our
// calling thresholds, an so don't appear in the Variant's alternate_bases()
// list. Such reads are added to the supporting read map keyed to this string
// value to indicate that they don't support reference but don't support an
// alternate allele either.
extern const char* const kSupportingUncalledAllele;

// Constants for the AD (depth by allele), DP (total depth), VAF (variant
// allele fraction), MF (methylation fraction), and MD (methylation depth)
// format fields.
extern const char* const kDPFormatField;
extern const char* const kADFormatField;
extern const char* const kVAFFormatField;
extern const char* const kMFFormatField;
extern const char* const kMDFormatField;

// Implements the less functionality needed to use an Allele as a key in a map.
struct OrderAllele {
  bool operator()(const Allele& allele1, const Allele& allele2) const {
    // Note we ignore count (and other potential fields) because they aren't
    // relevant in uses of this map.
    if (allele1.type() != allele2.type()) {
      return allele1.type() < allele2.type();
    } else {
      return allele1.bases() < allele2.bases();
    }
  }
};
using AlleleMap = std::map<Allele, std::string, OrderAllele>;

// Helper struct to store allele and position.
struct AlleleAtPosition {
  std::string alt_bases;
  AlleleType type;
  int position;
};

// Input options for SelectAltAlleles() function.
struct SelectAltAllelesInputOptions {
  const absl::node_hash_map<std::string, AlleleCount>&
      allele_counts_by_sample;
  bool create_complex_alleles;
  int prev_deletion_end;  // This is the end position of the previous deletion.
  // It is used to skip complex variant creation for alleles that are
  // overlapped by the previous deletion.
};

// Output options for SelectAltAlleles() function.
struct SelectAltAllelesResult {
  std::vector<Allele> alt_alleles;
  absl::node_hash_map<std::string, AlleleCount> allele_counts_mod;
  bool complex_variant_created;
  std::string ref_bases;
};

// Input options for SelectAltAllelesWithComplexVariant() function.
struct SelectAltAllelesWithComplexVariantInputOptions {
  const absl::node_hash_map<std::string,
                            AlleleCount>& allele_counts_by_sample;
  const std::vector<Allele>& alt_alleles;
};


// Helper function to generate a map of complex alleles to supporting reads.
// This function is used by SelectAltAllelesWithComplexVariant and it is
// defined in the header so that it can be tested.
absl::flat_hash_map<std::string, std::vector<std::string>>
CreateComplexAllelesSupport(
    const absl::flat_hash_map<std::string, std::vector<AlleleAtPosition>>&
        read_to_alt_alleles,
    int del_start, int del_len, absl::string_view del_allele_ref_bases);

// A very simple but highly sensitive variant caller.
//
// This class implements a very simple variant caller using the data
// in an AlleleCount proto.  It considers the distribution of Alleles
// observed at a position, and fills in the variant field with the Variant
// proto if there's reasonable evidence for a variant being at that site.
//
// The evidence standard is pretty loose: any allele that has at least
// min_count occurrences and that count is at least min_fraction
// of the total allele count at the site will be called.
//
// The heavy-lifting part of this code is just getting the alleles correct,
// in the case where there are multiple candidate alleles and the observed
// alleles need to be converted to their minimal VCF representation. For more
// information, see:
//
//   https://samtools.github.io/hts-specs/VCFv4.2.pdf
//
// Note that when multiple alleles satisfy our requirements, a multi-allelic
// Variant proto will be emitted.
//
// If a variant is emitted, the Variant proto created has only the minimal
// information needed to describe the call:
//
//  reference_name
//  start
//  end
//  reference_bases and alternate_bases
//
// No genotyping is attempted, so the variant.calls field is empty
// in the emitted Variant protos.
class VariantCaller {
 public:
  explicit VariantCaller(const VariantCallerOptions& options)
      : options_(options),
        sampler_(options.fraction_reference_sites_to_emit(),
                 options.random_seed()) {
    CHECK_GE(options_.min_count_snps(), 0) << "min_count_snps must be >= 0";
    CHECK_GE(options_.min_count_indels(), 0) << "min_count_indels must be >= 0";
    CHECK_GE(options_.min_fraction_snps(), 0.0)
        << "min_fraction_snps must be >= 0.0";
    CHECK_GE(options_.min_fraction_indels(), 0.0)
        << "min_fraction_indels must be >= 0.0";
    CHECK_GE(options_.fraction_reference_sites_to_emit(), 0.0)
        << "fraction_reference_sites_to_emit must be >= 0.0";
    CHECK_GE(options_.p_error(), 0.0) << "p_error must be >= 0.0";
    CHECK_GE(options_.max_gq(), 0) << "max_gq must be >= 0";
    CHECK_GE(options_.gq_resolution(), 0) << "gq_resolution must be >= 0";
    CHECK_GE(options_.ploidy(), 1) << "ploidy must be >= 1";
  }

  // MyClass is neither copyable nor movable.
  VariantCaller(const VariantCaller&) = delete;
  VariantCaller& operator=(const VariantCaller&) = delete;

  // Implements the filtering logic for a single allele. Multiple thresholds are
  // applied to determine if an allele is good enough to call. If the allele
  // passes all the thresholds, the function returns true. Otherwise, it is
  // rejected and the function returns false.
  bool AlleleFilter(const Allele& allele,
                    const AlleleCount& target_sample_allele_count,
                    absl::Span<const AlleleCount> all_samples_allele_counts,
                    absl::Span<const Allele> non_target_sample_alleles) const;
  // High-level API for calling variants in a region.
  //
  // Generate DeepVariantCall candidates for each position of the window.
  // AlleleCount objects from all samples are processed together per position.
  // Candidates are generated for those positions where there is enough support
  // for a candidate. There are two steps:
  // * Candidate generation attempt is made for the target sample.
  // * If candidate could not be generated in the first step due to not enough
  //   read support then another attempt is made to generate candidate from all
  //   the reads of all the samples.
  // Logic is implemented in AlleleFilter() function.
  std::vector<DeepVariantCall> CallsFromAlleleCounts(
      const std::unordered_map<std::string, AlleleCounter*>&
          allele_counts_wrapper,
      const std::string& target_sample);

  // High-level API for calculating potential variant position in a region.
  // This function is almost identical to CallsFromAlleleCounts except it
  // only calculates candidate positions.
  std::vector<int> CallPositionsFromAlleleCounts(
      const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
      const std::string& target_sample);

  // Iterates allele_counts for all samples and calls specified function F for
  // each candidate. Currently there are 2 use case: generate candidates,
  // generate candidate positions.
  template <class T>
  std::vector<T> AlleleCountsGenerator(
      std::optional<T> (VariantCaller::*F)(
          const absl::node_hash_map<std::string, AlleleCount>&,
          std::vector<AlleleCount>::const_iterator*,
          int& skip_next_count,
          int& prev_deletion_end) const) const;
  // Primary interface function for calling variants.
  //
  // Looks at the alleles in the provided AlleleCount proto and returns
  // either properly-formatted Variant proto specifying the non-reference call
  // or nullopt, indicating that the AlleleCount didn't meet the criteria for a
  // non-reference call. If a Variant is returned, it will have reference_name,
  // start, end set according to AlleleCount's position, reference_bases and
  // alternate_bases set based on the alleles in allele_count, along with an
  // appropriate end. The genotypes of the VariantCall will be set to -1 and -1
  // (diploid no-call).
  std::optional<DeepVariantCall> CallVariant(
      const absl::node_hash_map<std::string, AlleleCount>&
          allele_counts_per_sample,
      std::vector<AlleleCount>::const_iterator*
          target_sample_allele_count_iterator,
    int& skip_next_count,
    int& prev_deletion_end) const;

  // Adds supporting reads to the DeepVariantCall.
  void AddSupportingReads(
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      const AlleleMap& allele_map, absl::string_view target_sample,
      DeepVariantCall* call) const;

  // Adds allele counts in window around the position of the DeepVariantCall.
  void AddAdjacentAlleleFractionsAtPosition(
      int window_size,
      std::vector<AlleleCount>::const_iterator
          target_sample_allele_count_iterator,
      DeepVariantCall* call) const;

  // Helper function to combine the methylated reference sites and keep only the
  // positive strands.
  // For 5mC methylation, Pacbio marks only forward positions,
  // and the reverse is assumed.
  void MergeMethylatedAlleleCounts(
      const std::unordered_map<std::string, AlleleCounter*>& allele_counters)
      const;

  // Helper function to check if the allele is a reference site.
  // for methylation-aware phasing.
  bool IsReferenceSite(const AlleleCount& allele_count) const;

  // Helper function to check if the chromosome should be excluded from
  // methylation-aware phasing.
  bool IsExcludedMethylationContig(const std::string& chrom) const;

  // Computes and adds methylation fraction (MF) and methylation depth (MD)
  // to the INFO field of the variant
  void ComputeMethylationStats(
      const AlleleCount& target_sample_allele_count,
      const AlleleMap& allele_map, Variant* variant) const;

  void Clear() {
    for (auto& [sample, allele_counter] : allele_counters_per_sample_) {
      if (allele_counter != nullptr) {
        delete allele_counter;
      }
    }
  }

  // Helper function to create a VariantCaller for testing.
  static std::unique_ptr<VariantCaller> MakeTestVariantCallerFromAlleleCounts(
      VariantCallerOptions options, absl::Span<const AlleleCount> allele_counts,
      const std::string& sample_name) {
    auto caller = std::make_unique<VariantCaller>(options);
    caller->target_sample_ = sample_name;
    AlleleCounter* allele_counter =
        AlleleCounter::InitFromAlleleCounts(allele_counts);
    caller->allele_counters_per_sample_[sample_name] = allele_counter;
    return caller;
  }

 private:
  enum AlleleRejectionAcceptance {
    ACCEPTED,
    REJECTED_REF,
    REJECTED_LOW_SUPPORT,
    REJECTED_LOW_RATIO,
    REJECTED_OTHER
  };

  int min_count(const Allele& allele) const {
    return allele.type() == AlleleType::SUBSTITUTION
               ? options_.min_count_snps()
               : options_.min_count_indels();
  }
  double min_fraction(const Allele& allele) const {
    return allele.type() == AlleleType::SUBSTITUTION
               ? options_.min_fraction_snps()
               : options_.min_fraction_indels();
  }

  SelectAltAllelesResult SelectAltAlleles(
      const SelectAltAllelesInputOptions& options
  ) const;

  // Generate complex variants if there are multiple alleles and one of them is
  // a deletion and other alleles are not deletions. In all other cases return
  // SelectAltAlleles(). If complex alleles are generated then allele counts are
  // modified to reflect the new alleles.
  // It is expected that allele_counters_per_sample_ map contains the
  // target_sample_ key. The function will fail if it is not the case.
  SelectAltAllelesResult SelectAltAllelesWithComplexVariant(
      const SelectAltAllelesWithComplexVariantInputOptions& options
  ) const;

  AlleleRejectionAcceptance IsGoodAltAlleleWithReason(
      const Allele& allele, int total_count, bool apply_trio_coefficient) const;
  bool KeepReferenceSite() const;

  // This function duplicates functionality of CallVariant() to determine if
  // a position contains a candidate. If candidate conditions are met then
  // function returns a position of the candidate.
  std::optional<int> CallVariantPosition(
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      std::vector<AlleleCount>::const_iterator*
          target_sample_allele_count_iterator,
      int& skip_next_count,
      int& prev_deletion_end) const;

  const VariantCallerOptions options_;

  // Fraction of non-variant sites to emit as DeepVariantCalls.
  mutable nucleus::FractionalSampler sampler_;

  // AllelCounters hold a large internal vector, so we want to avoid copying
  // them.
  std::unordered_map<std::string, AlleleCounter*> allele_counters_per_sample_;
  // Name of the target sample for multi-allelic variant calling.
  std::string target_sample_;

  // Helper functions for methylation-aware phasing methylation merging
  // (PacBio only).
  // Returns true if the previous site is immediately before the current site.
  bool IsAdjacent(const AlleleCount& prev, const AlleleCount& curr) const;

  // Returns true if the site has a C reference or any C alternate allele.
  bool HasCRefOrAlt(const AlleleCount& allele_count) const;

  // Extracts and clears methylation from G alleles. Returns (read_key, level)
  // pairs.
  std::vector<std::tuple<std::string, int32_t, bool>>
  ExtractAndClearGSiteMethylation(AlleleCount& g_site) const;

  // Transfers methylation to the corresponding read’s C allele at the
  // previous site.
  void TransferMethylationToPrevC(
      AlleleCount& prev_allele_count,
      const std::vector<std::tuple<std::string, int32_t, bool>>&
          methylated_reads) const;

  FRIEND_TEST(VariantCallingTest, TestMultiAlleleWithDeletion);
  FRIEND_TEST(ComplexVariantTest, ComplexVariantTestCases);
  FRIEND_TEST(VariantCallingTest,
     TestCallVariantAddAdjacentAlleleFractionsAtPositionSize5);
  FRIEND_TEST(VariantCallingTest,
     TestCallVariantAddAdjacentAlleleFractionsAtPositionSize3);
  FRIEND_TEST(VariantCallingTest,
              TestCallVariantAddAdjacentAlleleFractionsAtPositionSize0);
  FRIEND_TEST(VariantCallingTest, TestRefSitesFraction);
  FRIEND_TEST(VariantCallingTest, TestCallVariantNew);
  friend class VariantCallingTest;
};

// Helper function
// If there are multiple deletions with different anchors at the same location
// this functions determines the deletions with the highest reads support and
// deletes all other deletions from the allele_map. In all other cases
// allele_map is not modified.
AlleleMap RemoveInvalidDels(const AlleleMap& allele_map,
                            absl::string_view ref_bases);

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_MULTISAMPLE_H_
