/*
 * Copyright 2017 Google Inc.
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
#ifndef LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_H_
#define LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_H_

#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/samplers.h"
#include "tensorflow/core/lib/gtl/optional.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// The alternate allele string for the gVCF "any" alternate allele.
extern const char *const kGVCFAltAllele;

// In a DeepVariantCall, reads can support an allele that didn't pass our
// calling thresholds, an so don't appear in the Variant's alternate_bases()
// list. Such reads are added to the supporting read map keyed to this string
// value to indicate that they don't support reference but don't support an
// alternate allele either.
extern const char* const kSupportingUncalledAllele;

// Constants for the AD (depth by allele), DP (total depth), and VAF (variant
// allele fraction) format fields.
extern const char *const kDPFormatField;
extern const char *const kADFormatField;
extern const char *const kVAFFormatField;

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
  // These functions invokes the CallVariant methods on each AlleleCount in
  // either the AlleleCounter object itself (via a call to Counts()) or on a
  // vector of AlleleCounts directly. This code processes each AlleleCount in
  // order, collecting up the DeepVariantCall protos at each site that
  // CallVariant says is a candidate variant. These DeepVariantCall protos are
  // returned in order.
  std::vector<DeepVariantCall> CallsFromAlleleCounter(
      const AlleleCounter& allele_counter) const;
  std::vector<DeepVariantCall> CallsFromAlleleCounts(
    const std::vector<AlleleCount>& allele_counts) const;

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
  tensorflow::gtl::optional<DeepVariantCall> CallVariant(
      const AlleleCount& allele_count) const;

 private:
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

  std::vector<Allele> SelectAltAlleles(const AlleleCount& allele_count) const;
  bool IsGoodAltAllele(const Allele& allele, const int total_count) const;
  bool KeepReferenceSite() const;

  const VariantCallerOptions options_;

  // Fraction of non-variant sites to emit as DeepVariantCalls.
  mutable nucleus::PhiloxFractionalSampler sampler_;
};


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_H_
