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

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include "absl/container/btree_map.h"
#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/statusor.h"
#include "absl/log/log.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace vcf_candidate_importer {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using tensorflow::string;

// Declared in .h.
const char* const kGVCFAltAllele = "<*>";
const char* const kSupportingUncalledAllele = "UNCALLED_ALLELE";
const char* const kDPFormatField = "DP";
const char* const kADFormatField = "AD";
const char* const kVAFFormatField = "VAF";

// The VCF/Variant allele string to use when you don't have any alt alleles.
const char* const kNoAltAllele = ".";

namespace {

template <class T>
std::vector<T> AsVector(const google::protobuf::RepeatedPtrField<T>& container) {
  return std::vector<T>(container.begin(), container.end());
}

// Adds a single VariantCall with sample_name, genotypes, and gq (bound to the
// "GQ" key of info with a numerical value of gq, if provided) to variant.
void AddGenotypes(const string& sample_name, const std::vector<int>& genotypes,
                  Variant* variant) {
  CHECK(variant != nullptr);

  VariantCall* call = variant->add_calls();
  call->set_call_set_name(sample_name);
  for (const auto genotype : genotypes) {
    call->add_genotype(genotype);
  }
}

void FillVariant(const string& reference_name, int variant_start,
                 const string& ref_bases, const string& sample_name,
                 const std::vector<std::string>& alternate_bases,
                 Variant* variant) {
  variant->set_reference_name(reference_name);
  variant->set_start(variant_start);
  variant->set_reference_bases(ref_bases);
  variant->set_end(variant->start() + ref_bases.size());
  AddGenotypes(sample_name, {-1, -1}, variant);

  for (const string& alt : alternate_bases) {
    variant->add_alternate_bases(alt);
  }
}

void MakeVariantConsistentWithRefAndAlts(const string& refbases,
                                         const std::vector<Allele>& alt_alleles,
                                         Variant* variant_to_fix) {
  if (variant_to_fix->reference_bases() == refbases) {
    // No fix needed if the reference bases are identical.
    return;
  }
  QCHECK_NE(variant_to_fix->reference_bases().length(), refbases.length())
      << "Proposed variant has incorrect ref bases: "
      << "Problematic variant=" << variant_to_fix->DebugString();

  if (variant_to_fix->reference_bases().length() < refbases.length()) {
    QCHECK(absl::StartsWith(refbases, variant_to_fix->reference_bases()))
        << "Proposed variant has incorrect ref bases: "
        << "Problematic variant=" << variant_to_fix->DebugString();
    string suffix = refbases.substr(variant_to_fix->reference_bases().length(),
                                    refbases.length());
    variant_to_fix->set_reference_bases(
        absl::StrCat(variant_to_fix->reference_bases(), suffix));
    for (int i = 0; i < variant_to_fix->alternate_bases_size(); ++i) {
      variant_to_fix->set_alternate_bases(
          i, absl::StrCat(variant_to_fix->alternate_bases(i), suffix));
    }
    variant_to_fix->set_end(variant_to_fix->end() + suffix.length());
  }
}

// Assumption: `short_str` is a prefix of `long_str`.
// Return the suffix on long_str.
string GetSuffixFromTwoAlleles(const string& short_str,
                               const string& long_str) {
  QCHECK(absl::StartsWith(long_str, short_str))
      << short_str << " has to be a prefix of " << long_str;
  return long_str.substr(short_str.length(), long_str.length());
}

}  // namespace

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
string CalcRefBases(const string& ref_bases,
                    const std::vector<Allele>& alt_alleles) {
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
    return absl::StrCat(ref_bases, max_elt->bases().substr(1));
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
string MakeAltAllele(const string& prefix, const string& variant_ref,
                     const uint32_t from) {
  const auto postfix =
      from >= variant_ref.length() ? "" : variant_ref.substr(from);
  return absl::StrCat(prefix, postfix);
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
                   const string& allele_map_refbases, Variant* variant) {
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

    absl::btree_map<std::string, const Allele*, std::less<>> alt_to_alleles;
    for (const auto& entry : allele_map) {
      const string key = SimplifyRefAlt(allele_map_refbases, entry.second);
      alt_to_alleles[key] = &entry.first;
    }
    CHECK(alt_to_alleles.size() == allele_map.size())
        << "Non-unique alternative alleles!";
    for (const string& alt : variant->alternate_bases()) {
      const string simplified_ref_alt =
          SimplifyRefAlt(variant->reference_bases(), alt);
      int count_of_allele = 0;
      auto found = alt_to_alleles.find(simplified_ref_alt);
      if (found != alt_to_alleles.end()) {
        count_of_allele = (*found->second).count();
      }
      double this_vaf = 0.0;
      if (dp > 0) {
        this_vaf = 1.0 * count_of_allele / dp;
      }
      ad.push_back(count_of_allele);
      vaf.push_back(this_vaf);
    }

    nucleus::SetInfoField(kADFormatField, ad, call);
    nucleus::SetInfoField(kVAFFormatField, vaf, call);
  }
}

// Returns true if the current site should be emited, even if it's a reference
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
    std::optional<DeepVariantCall> call = CallVariant(allele_count);
    if (call) {
      variants.push_back(*call);
    }
  }

  return variants;
}

bool is_uncalled_genotype(const Variant& variant) {
  if (variant.calls_size() >= 1) {
    VariantCall call = variant.calls().Get(0);
    if (call.genotype().size() >= 2) {
      return call.genotype().Get(0) == -1 && call.genotype().Get(1) == -1;
    }
  }
  return false;
}

std::vector<DeepVariantCall> VariantCaller::CallsFromVcf(
    const std::vector<AlleleCount>& allele_counts,
    const Range& range,
    nucleus::VcfReader* vcf_reader_ptr) const {
  std::vector<Variant> variants_in_region;
  nucleus::StatusOr<std::shared_ptr<nucleus::VariantIterable>> status =
      vcf_reader_ptr->Query(range);
  if (status.ok()) {
    std::shared_ptr<nucleus::VariantIterable> variants = status.ValueOrDie();
    bool warn_missing = false;
    for (const auto& v : variants) {
      const Variant* variant = v.ValueOrDie();
      // This ensures we only keep variants that start in this region.
      // By default, vcf_reader->Query() returns all variants that overlap a
      // region, which can incorrectly cause the same variant to be processed
      // multiple times.
      if (variant->start() >= range.start()) {
        if (options_.skip_uncalled_genotypes() &&
            is_uncalled_genotype(*variant)) {
          if (!warn_missing) {
            LOG(WARNING) << "Uncalled genotypes (./.) present in VCF. These "
                            "are skipped.";
            warn_missing = true;
          }
          continue;
        }
        Variant clean_variant;
        FillVariant(variant->reference_name(), variant->start(),
                    variant->reference_bases(), options_.sample_name(),
                    AsVector<std::string>(variant->alternate_bases()),
                    &clean_variant);
        variants_in_region.push_back(clean_variant);
      }
    }
  } else if (status.error_message() == "Cannot query without an index") {
    LOG(FATAL) << "Error in VariantCaller::CallsFromVcf: "
               << status.error_message();
  } else {
    LOG(WARNING)
        << nucleus::MakeIntervalStr(range)
        << " cannot be found in proposed VCF header. Skip this region.";
  }
  return CallsFromVariantsInRegion(allele_counts, variants_in_region);
}

std::vector<int> VariantCaller::CallPositionsFromVcf(
    const std::vector<AlleleCount>& allele_counts, const Range& range,
    nucleus::VcfReader* vcf_reader_ptr) const {
  std::vector<Variant> variants_in_region;
  std::vector<int> positions;
  nucleus::StatusOr<std::shared_ptr<nucleus::VariantIterable>> status =
      vcf_reader_ptr->Query(range);
  if (status.ok()) {
    std::shared_ptr<nucleus::VariantIterable> variants = status.ValueOrDie();
    bool warn_missing = false;
    for (const auto& v : variants) {
      const Variant* variant = v.ValueOrDie();
      // This ensures we only keep variants that start in this region.
      // By default, vcf_reader->Query() returns all variants that overlap a
      // region, which can incorrectly cause the same variant to be processed
      // multiple times.
      if (variant->start() >= range.start()) {
        if (options_.skip_uncalled_genotypes() &&
            is_uncalled_genotype(*variant)) {
          if (!warn_missing) {
            LOG(WARNING) << "Uncalled genotypes (./.) present in VCF. These "
                            "are skipped.";
            warn_missing = true;
          }
          continue;
        }
        // This is a good variant, save the position.
        positions.push_back(variant->start());
      }
    }
  } else if (status.error_message() == "Cannot query without an index") {
    LOG(FATAL) << "Error in VariantCaller::CallsFromVcf: "
               << status.error_message();
  } else {
    LOG(WARNING)
        << nucleus::MakeIntervalStr(range)
        << " cannot be found in proposed VCF header. Skip this region.";
  }
  return positions;
}

std::vector<DeepVariantCall> VariantCaller::CallsFromVcf(
    const AlleleCounter& allele_counter,
    nucleus::VcfReader* vcf_reader_ptr) const {
  return CallsFromVcf(allele_counter.Counts(), allele_counter.Interval(),
                      vcf_reader_ptr);
}

std::vector<int> VariantCaller::CallPositionsFromVcf(
    const AlleleCounter& allele_counter,
    nucleus::VcfReader* vcf_reader_ptr) const {
  return CallPositionsFromVcf(allele_counter.Counts(),
                              allele_counter.Interval(),
                              vcf_reader_ptr);
}

std::vector<DeepVariantCall> VariantCaller::CallsFromVariantsInRegion(
    const std::vector<AlleleCount>& allele_counts,
    const std::vector<Variant>& variants_in_region) const {
  std::vector<DeepVariantCall> calls;
  // For each variant in the region, loop through AlleleCounts to find a match
  // to the variant position. At each match, add the supporting reads.
  for (const auto& v : variants_in_region) {
    std::optional<DeepVariantCall> call = ComputeVariant(v, allele_counts);
    if (call) {
      calls.push_back(*call);
    }
  }
  return calls;
}

std::optional<DeepVariantCall> VariantCaller::ComputeVariant(
    const Variant& variant,
    const std::vector<AlleleCount>& allele_counts) const {
  DeepVariantCall call;
  *call.mutable_variant() = variant;
  Variant* m_variant = call.mutable_variant();
  AlleleCount allele_count_match;

  int idx = AlleleIndex(allele_counts, variant.start());
  if (idx != -1) {
    allele_count_match = allele_counts[idx];
    if (!nucleus::AreCanonicalBases(allele_count_match.ref_base())) {
      // We don't emit calls at any site in the genome that isn't one of the
      // canonical DNA bases (one of A, C, G, or T).
      return std::nullopt;
    }
  }
  // If idx=-1 and no allele count matches we proceed with
  // an empty allele_count_match object which is used to help return
  // a missing genotype with no observed evidence.

  std::vector<Allele> alt_alleles = SelectAltAlleles(allele_count_match);
  string refbases = CalcRefBases(allele_count_match.ref_base(), alt_alleles);
  MakeVariantConsistentWithRefAndAlts(refbases, alt_alleles, m_variant);

  // Compute the map from read alleles to the alleles we'll use in our Variant.
  // Add the alternate alleles from our allele_map to the variant.
  const AlleleMap allele_map =
      BuildAlleleMap(allele_count_match, alt_alleles, refbases);

  AddReadDepths(allele_count_match, allele_map, refbases, m_variant);
  AddSupportingReads(allele_count_match.read_alleles(), allele_map, refbases,
                     &call);
  return std::make_optional(call);
}

std::optional<DeepVariantCall> VariantCaller::CallVariant(
    const AlleleCount& allele_count) const {
  if (!nucleus::AreCanonicalBases(allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return std::nullopt;
  }

  std::vector<Allele> alt_alleles = SelectAltAlleles(allele_count);
  if (alt_alleles.empty() && !KeepReferenceSite()) {
    return std::nullopt;
  }
  const string refbases = CalcRefBases(allele_count.ref_base(), alt_alleles);
  std::vector<std::string> alternate_bases;
  // Compute the map from read alleles to the alleles we'll use in our Variant.
  // Add the alternate alleles from our allele_map to the variant.
  const AlleleMap allele_map =
      BuildAlleleMap(allele_count, alt_alleles, refbases);
  for (const auto& elt : allele_map) {
    alternate_bases.push_back(elt.second);
  }
  // If we don't have any alt_alleles, we are generating a reference site so
  // add in the kNoAltAllele.
  if (alt_alleles.empty()) alternate_bases.push_back(kNoAltAllele);
  std::sort(alternate_bases.begin(), alternate_bases.end());

  DeepVariantCall call;
  Variant* variant = call.mutable_variant();
  string sample_name = options_.sample_name();
  if (variant->calls_size() > 0 && !variant->calls(0).call_set_name().empty()) {
    sample_name = variant->calls(0).call_set_name();
  }
  // Creates a non-reference Variant proto based on the information in
  // allele_count and alt_alleles. This variant starts at the position of
  // allele_count with the same reference_name. The reference_bases are
  // calculated based on the alt_alleles, which are also set appropriately for
  // the variant. For convenience, the alt_alleles are sorted.
  FillVariant(allele_count.position().reference_name(),
              allele_count.position().position(), refbases, sample_name,
              alternate_bases, variant);
  AddReadDepths(allele_count, allele_map, refbases, variant);
  AddSupportingReads(allele_count.read_alleles(), allele_map, refbases, &call);
  return std::make_optional(call);
}

void VariantCaller::AddSupportingReads(
    const ::google::protobuf::Map<std::string, Allele>& read_alleles,
    const AlleleMap& allele_map, const string& refbases,
    DeepVariantCall* call) const {
  string suffix = "";
  if (call->variant().reference_bases().length() > refbases.length()) {
    suffix =
        GetSuffixFromTwoAlleles(refbases, call->variant().reference_bases());
  }
  // Iterate over each read in the allele_count, and add its name to the
  // supporting reads of for the Variant allele it supports.
  const string unknown_allele = kSupportingUncalledAllele;
  for (const auto& read_name_allele : read_alleles) {
    const string& read_name = read_name_allele.first;
    const Allele& allele = read_name_allele.second;

    // Skip reference supporting reads, as they aren't included in the
    // supporting reads for alternate alleles.
    if (allele.type() != AlleleType::REFERENCE) {
      auto it = allele_map.find(allele);
      const string supported_allele = it == allele_map.end()
                                          ? unknown_allele
                                          : absl::StrCat(it->second, suffix);
      DeepVariantCall_SupportingReads& supports =
          (*call->mutable_allele_support())[supported_allele];
      supports.add_read_names(read_name);
      DeepVariantCall_SupportingReadsExt& support_infos =
          (*call->mutable_allele_support_ext())[supported_allele];
      DeepVariantCall_ReadSupport* read_info = support_infos.add_read_infos();
      read_info->set_read_name(read_name);
      read_info->set_is_low_quality(allele.is_low_quality());
    } else if (options_.track_ref_reads()) {
      call->add_ref_support(read_name);
      DeepVariantCall_SupportingReadsExt& support_infos =
          (*call->mutable_ref_support_ext());
      DeepVariantCall_ReadSupport* read_info = support_infos.add_read_infos();
      read_info->set_read_name(read_name);
      read_info->set_is_low_quality(allele.is_low_quality());
    }
  }
}

}  // namespace vcf_candidate_importer
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
