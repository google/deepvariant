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

#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/node_hash_map.h"
#include "absl/strings/string_view.h"
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
using learning::genomics::deepvariant::AlleleType;
using learning::genomics::deepvariant::DeepVariantCall;
using learning::genomics::deepvariant::VariantCallerOptions;

// The alternate allele string for the gVCF "any" alternate allele.
extern const char* const kGVCFAltAllele;

// In a DeepVariantCall, reads can support an allele that didn't pass our
// calling thresholds, an so don't appear in the Variant's alternate_bases()
// list. Such reads are added to the supporting read map keyed to this string
// value to indicate that they don't support reference but don't support an
// alternate allele either.
extern const char* const kSupportingUncalledAllele;

// Constants for the AD (depth by allele), DP (total depth), and VAF (variant
// allele fraction) format fields.
extern const char* const kDPFormatField;
extern const char* const kADFormatField;
extern const char* const kVAFFormatField;

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
  // Logic is implemented in SelectAltAlleles() function.
  std::vector<DeepVariantCall> CallsFromAlleleCounts(
      const std::unordered_map<std::string, AlleleCounter*>&
          allele_counts_wrapper,
      const std::string& target_sample) const;

  // High-level API for calculating potential variant position in a region.
  // This function is almost identical to CallsFromAlleleCounts except it
  // only calculates candidate positions.
  std::vector<int> CallPositionsFromAlleleCounts(
      const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
      const std::string& target_sample) const;

  // Iterates allele_counts for all samples and calls specified function F for
  // each candidate. Currently there are 2 use case: generate candidates,
  // generate candidate positions.
  template <class T>
  std::vector<T> AlleleCountsGenerator(
      const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
      const std::string& target_sample,
      std::optional<T> (VariantCaller::*F)(
          const absl::node_hash_map<std::string, AlleleCount>&,
          const std::string&, const std::vector<AlleleCount>*,
          std::vector<AlleleCount>::const_iterator*) const) const;
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
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      const std::string& target_sample,
      const std::vector<AlleleCount>* target_sample_allele_counts = nullptr,
      std::vector<AlleleCount>::const_iterator*
          target_sample_allele_count_iterator = nullptr) const;

  // Adds supporting reads to the DeepVariantCall.
  void AddSupportingReads(
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      const AlleleMap& allele_map, const std::string& target_sample,
      DeepVariantCall* call) const;

  // Adds allele counts in window around the position of the DeepVariantCall.
  void AddAdjacentAlleleFractionsAtPosition(
      int window_size,
      const std::vector<AlleleCount>& target_sample_allele_counts,
      std::vector<AlleleCount>::const_iterator
          target_sample_allele_count_iterator,
      DeepVariantCall* call) const;

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

  std::vector<Allele> SelectAltAlleles(
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      absl::string_view target_sample) const;
  AlleleRejectionAcceptance IsGoodAltAlleleWithReason(
      const Allele& allele, const int total_count,
      const bool apply_trio_coefficient) const;
  bool KeepReferenceSite() const;

  // This function duplicates functionality of CallVariant() to determine if
  // a position contains a candidate. If candidate conditions are met then
  // function returns a position of the candidate.
  std::optional<int> CallVariantPosition(
      const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
      const std::string& target_sample,
      const std::vector<AlleleCount>* target_sample_allele_counts = nullptr,
      std::vector<AlleleCount>::const_iterator*
          target_sample_allele_count_iterator = nullptr) const;

  const VariantCallerOptions options_;

  // Fraction of non-variant sites to emit as DeepVariantCalls.
  mutable nucleus::FractionalSampler sampler_;
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
