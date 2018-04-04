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

#include "deepvariant/variant_calling.h"

#include <algorithm>
#include <numeric>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

namespace {
// Used for sorting RepeatedPtrField below.
struct StringPtrLessThan {
    bool operator() (const string* x, const string* y) const {
          return *x < *y;
    }
};
}

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using tensorflow::gtl::nullopt;
using tensorflow::gtl::optional;
using tensorflow::gtl::make_optional;
using tensorflow::strings::StrCat;

// Declared in .h.
const char* const kGVCFAltAllele = "<*>";
const char* const kSupportingUncalledAllele = "UNCALLED_ALLELE";
const char *const kDPFormatField = "DP";
const char *const kADFormatField = "AD";
const char *const kVAFFormatField = "VAF";

// The VCF/Variant allele string to use when you don't have any alt alleles.
const char* const kNoAltAllele = ".";

////////////////////////////////////////////////////////////////////////////////
//
// Code for creating non-reference Variant calls.
//
////////////////////////////////////////////////////////////////////////////////

// Get the 'deletion' size of allele, which is the length of the
// bases if allele is a deletion, or -1 otherwise.  A helper
// function for CalcRefBases.
int DeletionSize(const Allele& allele) {
  return allele.type() == AlleleType::DELETION ? allele.bases().length() : -1;
}

// Get the bases to use as the reference bases in a Variant proto.
//
// The reference bases in a variant proto represent the longest substitution
// of bases on the reference genome needed to describe a substitution by
// one of alt_alleles in a sample. What this means is that if alt_alleles
// doesn't include any deletions, this is simply the reference bases of our
// AlleleCount. But if one of the alt_alleles is a deletion, we need to
// use those bases as our reference.  And if there are multiple deletions
// at a site, we need to use the longest deletion allele.
string CalcRefBases(const string &ref_bases,
                    const std::vector<Allele> &alt_alleles) {
  if (alt_alleles.empty()) {
    // We don't have any alternate alleles, so used the provided ref_bases.
    return ref_bases;
  }

  const auto max_elt =
      std::max_element(alt_alleles.cbegin(), alt_alleles.cend(),
                       [](const Allele& allele1, const Allele& allele2) {
                         return DeletionSize(allele1) < DeletionSize(allele2);
                       });
  if (max_elt->type() != AlleleType::DELETION) {
    return ref_bases;
  } else {
    // Deletion alleles may have an anchor base that is the reference or some
    // other base, but a Variant must have a reference sequence that starts with
    // the reference base. The index 1 skips the first base of the deletion,
    // which is the anchor base of the deletion.
    CHECK(max_elt->bases().size() > 1)
        << "Saw invalid deletion allele with too few bases"
        << max_elt->ShortDebugString();
    return StrCat(ref_bases, max_elt->bases().substr(1));
  }
}

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
string MakeAltAllele(const string& prefix,
                     const string& variant_ref,
                     const uint32_t from) {
  const auto postfix = from >= variant_ref.length()
                       ? "" : variant_ref.substr(from);
  return StrCat(prefix, postfix);
}

// Is allele a good alternative allele for a Variant proto?
//
// A good alt allele is one that is a substitution, insertion, or deletion,
// and satisfies our min count and min fraction requirements.
bool VariantCaller::IsGoodAltAllele(const Allele& allele,
                                    const int total_count) const {
  return allele.type() != AlleleType::REFERENCE &&
         allele.type() != AlleleType::SOFT_CLIP &&
         allele.count() >= min_count(allele) &&
         (1.0 * allele.count()) / total_count >= min_fraction(allele);
}

// Select the subset of GoodAltAlleles from the alleles of allele_count.
//
// Returns the vector of allele objects from allele_count that satisfy
// IsGoodAltAllele().
std::vector<Allele> VariantCaller::SelectAltAlleles(
    const AlleleCount& allele_count) const {
  const std::vector<Allele> alleles = SumAlleleCounts(allele_count);
  const int total_count = TotalAlleleCounts(allele_count);
  std::vector<Allele> alt_alleles;
  for (const auto& allele : alleles) {
    if (IsGoodAltAllele(allele, total_count)) {
      alt_alleles.push_back(allele);
    }
  }
  return alt_alleles;
}

// Adds a single VariantCall with sample_name, genotypes, and gq (bound to the
// "GQ" key of info with a numerical value of gq, if provided) to variant.
void AddGenotypes(const string& sample_name,
                  const std::vector<int>& genotypes, Variant* variant) {
  CHECK(variant != nullptr);

  VariantCall* call = variant->add_calls();
  call->set_call_set_name(sample_name);
  for (const auto genotype : genotypes) {
    call->add_genotype(genotype);
  }
}

// Implements the less functionality needed to use an Allele as an key in a map.
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
using AlleleMap = std::map<Allele, string, OrderAllele>;

AlleleMap BuildAlleleMap(const AlleleCount& allele_count,
                         const std::vector<Allele>& alt_alleles,
                         const string& ref_bases) {
  AlleleMap allele_map;

  // Compute the alt alleles, recording the mapping from each Allele to its
  // corresponding allele in the Variant format.
  for (const auto& alt_allele : alt_alleles) {
    const string& alt_bases = alt_allele.bases();
    switch (alt_allele.type()) {
      case AlleleType::SUBSTITUTION:
      case AlleleType::INSERTION:
        allele_map[alt_allele] = MakeAltAllele(alt_bases, ref_bases, 1);
        break;
      case AlleleType::DELETION: {
        // The prefix base for a deletion should be the first base of the
        // deletion allele, which can be reference but might not be.
        CHECK(alt_bases.size() > 1)
            << "Saw invalid deletion allele with too few bases"
            << alt_allele.ShortDebugString();
        // The prefix base here is the anchor base of the deletion allele, which
        // is the first base of the alt_bases string.
        allele_map[alt_allele] =
            MakeAltAllele(alt_bases.substr(0, 1), ref_bases, alt_bases.size());
        break;
      }
      case AlleleType::SOFT_CLIP:
        // We don't want to add SOFT_CLIP alleles to our map.
        break;
      default:
        // this includes AlleleType::REFERENCE which should have been removed
        LOG(FATAL) << "Unexpected alt allele " << alt_allele.DebugString();
    }
  }

  return allele_map;
}

// Adds the DP, AD, and VAF VCF fields to the first VariantCall of Variant.
// DP: the total number of observed reads at the site.
// AD: the number of reads supporting each of our ref and alt alleles.
// VAF: the allele fraction of the variants (only including alt alleles).
// These are calculated from the provided allele_count information. The
// allele_map is needed to map between the Variant reference and alternate_bases
// and the Alleles used in allele_count.
void AddReadDepths(const AlleleCount& allele_count, const AlleleMap& allele_map,
                   Variant* variant) {
  // Set the DP to the total good reads seen at this position.
  VariantCall* call = variant->mutable_calls(0);
  nucleus::SetInfoField(kDPFormatField, TotalAlleleCounts(allele_count), call);

  if (variant->alternate_bases_size() == 1 &&
      (variant->alternate_bases(0) == kNoAltAllele ||
       variant->alternate_bases(0) == kGVCFAltAllele)) {
    // Variant has no alts or is a a gVCF record so only DP is meaningful.
    return;
  } else {
    int dp = TotalAlleleCounts(allele_count);
    // Build up AD and VAF.
    std::vector<int> ad;
    std::vector<double> vaf;
    ad.push_back(allele_count.ref_supporting_read_count());

    std::map<tensorflow::StringPiece, const Allele*> alt_to_alleles;
    for (const auto& entry : allele_map) {
      alt_to_alleles[entry.second] = &entry.first;
    }
    CHECK(alt_to_alleles.size() == allele_map.size())
        << "Non-unique alternative alleles!";
    for (const string& alt : variant->alternate_bases()) {
      const Allele& allele = *alt_to_alleles.find(alt)->second;
      ad.push_back(allele.count());
      vaf.push_back(1.0 * allele.count() / dp);
    }

    nucleus::SetInfoField(kADFormatField, ad, call);
    nucleus::SetInfoField(kVAFFormatField, vaf, call);
  }
}

// Returns true if the current site should be emited, even if its a reference
// site. This function is used to return reference site samples if the
// member variable fraction_reference_sites_to_emit >= 0.0 by pulling draws
// from a random number and returning true if the number <= the threshold.
bool VariantCaller::KeepReferenceSite() const {
  return options_.fraction_reference_sites_to_emit() > 0.0 && sampler_.Keep();
}

std::vector<DeepVariantCall> VariantCaller::CallsFromAlleleCounter(
    const AlleleCounter& allele_counter) const {
  return CallsFromAlleleCounts(allele_counter.Counts());
}

std::vector<DeepVariantCall> VariantCaller::CallsFromAlleleCounts(
    const std::vector<AlleleCount>& allele_counts) const {
  std::vector<DeepVariantCall> variants;
  for (const AlleleCount& allele_count : allele_counts) {
    optional<DeepVariantCall> call = CallVariant(allele_count);
    if (call) {
      variants.push_back(*call);
    }
  }

  return variants;
}

optional<DeepVariantCall> VariantCaller::CallVariant(
    const AlleleCount& allele_count) const {
  if (!nucleus::AreCanonicalBases(allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return nullopt;
  }

  const std::vector<Allele>& alt_alleles = SelectAltAlleles(allele_count);
  if (alt_alleles.empty() && !KeepReferenceSite()) {
    return nullopt;
  }
  // Creates a non-reference Variant proto based on the information in
  // allele_count and alt_alleles. This variant starts at the position of
  // allele_count with the same reference_name. The reference_bases are
  // calculated based on the alt_alleles, which are also set appropriately for
  // the variant. For convenience, the alt_alleles are sorted. Also adds a
  // single VariantCall to the Variant, with sample_name and uncalled diploid
  // genotypes.
  DeepVariantCall call;
  Variant* variant = call.mutable_variant();
  variant->set_reference_name(allele_count.position().reference_name());
  variant->set_start(allele_count.position().position());
  const string refbases = CalcRefBases(allele_count.ref_base(), alt_alleles);
  variant->set_reference_bases(refbases);
  variant->set_end(variant->start() + refbases.size());
  AddGenotypes(options_.sample_name(), {-1, -1}, variant);

  // Compute the map from read alleles to the alleles we'll use in our Variant.
  // Add the alternate alleles from our allele_map to the variant.
  AlleleMap allele_map = BuildAlleleMap(allele_count, alt_alleles, refbases);
  for (const auto& elt : allele_map) {
    variant->add_alternate_bases(elt.second);
  }
  // If we don't have any alt_alleles, we are generating a reference site so
  // add in the kNoAltAllele.
  if (alt_alleles.empty()) variant->add_alternate_bases(kNoAltAllele);
  std::sort(variant->mutable_alternate_bases()->pointer_begin(),
            variant->mutable_alternate_bases()->pointer_end(),
            StringPtrLessThan());

  AddReadDepths(allele_count, allele_map, variant);

  // Iterate over each read in the allele_count, and add its name to the
  // supporting reads of for the Variant allele it supports.
  const string unknown_allele = kSupportingUncalledAllele;
  for (const auto& read_name_allele : allele_count.read_alleles()) {
    const string& read_name = read_name_allele.first;
    const Allele& allele = read_name_allele.second;

    // Skip reference supporting reads, as they aren't included in the
    // supporting reads for alternate alleles.
    if (allele.type() != AlleleType::REFERENCE) {
      auto it = allele_map.find(allele);
      const string& supported_allele =
          it == allele_map.end() ? unknown_allele : it->second;
      auto& supports = (*call.mutable_allele_support())[supported_allele];
      supports.add_read_names(read_name);
    }
  }

  return make_optional(call);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
