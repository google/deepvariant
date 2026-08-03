/*
 * Copyright 2026 Google LLC.
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

// Shared utilities for variant calling used by both the single-sample
// (variant_calling.cc) and multi-sample (variant_calling_multisample.cc)
// variant callers.
//
#ifndef LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_UTILS_H_
#define LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_UTILS_H_

#include <cstdint>
#include <string>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/btree_map.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace variant_calling_utils {

// The alternate allele string for the gVCF "any" alternate allele.
extern const char* const kGVCFAltAllele;

// In a DeepVariantCall, reads can support an allele that didn't pass our
// calling thresholds, and so don't appear in the Variant's alternate_bases()
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

// The VCF/Variant allele string to use when you don't have any alt alleles.
extern const char* const kNoAltAllele;

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
using AlleleMap = absl::btree_map<Allele, std::string, OrderAllele>;

// Get the 'deletion' size of allele, which is the length of the
// bases if allele is a deletion, or -1 otherwise.  A helper
// function for CalcRefBases.
int DeletionSize(const Allele& allele);

// Get the bases to use as the reference bases in a Variant proto.
//
// The reference bases in a variant proto represent the longest substitution
// of bases on the reference genome needed to describe a substitution by
// one of alt_alleles in a sample. What this means is that if alt_alleles
// doesn't include any deletions, this is simply the reference bases of our
// AlleleCount. But if one of the alt_alleles is a deletion, we need to
// use those bases as our reference.  And if there are multiple deletions
// at a site, we need to use the longest deletion allele.
//
// If rejected_alleles is provided (used in multisample calling), it will also
// be considered when finding the longest deletion.
std::string CalcRefBases(absl::string_view ref_bases,
                         absl::Span<const Allele> alt_alleles,
                         absl::Span<const Allele> rejected_alleles = {});

// Constructs an alt allele from the prefix bases and the reference bases.
//
// This function helps create alt alleles for a variant proto. The complex logic
// here is to deal with the fact that the variant_ref bases aren't the simple
// single reference base context that the Allele objects are in but rather the
// actual reference bases of the variant, which could include a long series of
// bases if there's a deletion allele.
//
// This function takes a prefix of bases and concatenates those bases onto the
// appropriate substring of variant_ref. The substring starts at the from
// argument and runs to the end of variant_ref string, provided from isn't
// beyond the end of variant_ref.
//
// Suppose that we have variant_ref == "ACGT" due to a deletion, and our alleles
// are "C" [SNP] and "ATTT" [INSERTION] along with our "ACGT" [DELETION]. Each
// allele comes into this function with the following arguments:
//
//   "C" [SNP]    : prefix="C" and from=1
//   "ATTT" [INS] : prefix="ATTT" and from=1
//   "ACGT" [DEL] : prefix="A" (original ref base) and from=4
//
// This function will produce appropriate alleles that correct for the new
// reference bases due to the deletion as:
//
//   "C" [SNP]    => "C" + "CGT" => "CCGT", putting back deleted bases
//   "ATTT" [INS] => "ATTT" + "CGT" => "ATTTCGT", putting back deleted bases
//   "ACGT" [DEL] => "A" + "" (from >= "ACGT".length()) => "A"
//
std::string MakeAltAllele(absl::string_view prefix,
                          absl::string_view variant_ref, uint32_t from);

// Adds a single VariantCall with sample_name, genotypes, and gq (bound to the
// "GQ" key of info with a numerical value of gq, if provided) to variant.
void AddGenotypes(const std::string& sample_name,
                  absl::Span<const int> genotypes,
                  nucleus::genomics::v1::Variant* variant);

// Adds the DP, AD, and VAF VCF fields to the first VariantCall of Variant.
// DP: the total number of observed reads at the site.
// AD: the number of reads supporting each of our ref and alt alleles.
// VAF: the allele fraction of the variants (only including alt alleles).
// These are calculated from the provided allele_count information. The
// allele_map is needed to map between the Variant reference and alternate_bases
// and the Alleles used in allele_count.
void AddReadDepths(const AlleleCount& allele_count, const AlleleMap& allele_map,
                   absl::string_view allele_map_refbases,
                   nucleus::genomics::v1::Variant* variant);

}  // namespace variant_calling_utils
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_VARIANT_CALLING_UTILS_H_
