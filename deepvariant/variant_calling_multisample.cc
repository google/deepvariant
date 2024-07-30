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

#include "deepvariant/variant_calling_multisample.h"

#include <stdlib.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace multi_sample {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

// Declared in .h.
const char* const kGVCFAltAllele = "<*>";
const char* const kSupportingUncalledAllele = "UNCALLED_ALLELE";
const char* const kDPFormatField = "DP";
const char* const kADFormatField = "AD";
const char* const kVAFFormatField = "VAF";

// The VCF/Variant allele string to use when you don't have any alt alleles.
const char* const kNoAltAllele = ".";

namespace {
// Used for sorting RepeatedPtrField below.
struct StringPtrLessThan {
  bool operator()(const std::string* x, const std::string* y) const {
    return *x < *y;
  }
};
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
std::string CalcRefBases(absl::string_view ref_bases,
                         const std::vector<Allele>& alt_alleles) {
  if (alt_alleles.empty()) {
    // We don't have any alternate alleles, so used the provided ref_bases.
    return std::string(ref_bases);
  }

  const auto max_elt =
      std::max_element(alt_alleles.cbegin(), alt_alleles.cend(),
                       [](const Allele& allele1, const Allele& allele2) {
                         return DeletionSize(allele1) < DeletionSize(allele2);
                       });
  if (max_elt->type() != AlleleType::DELETION) {
    return std::string(ref_bases);
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
std::string MakeAltAllele(const std::string_view prefix,
                          const std::string& variant_ref, const uint32_t from) {
  const auto postfix =
      from >= variant_ref.length() ? "" : variant_ref.substr(from);
  return absl::StrCat(prefix, postfix);
}

// Is allele a good alternative allele for a Variant proto?
//
// A good alt allele is one that is a substitution, insertion, or deletion,
// and satisfies our min count and min fraction requirements.
VariantCaller::AlleleRejectionAcceptance
VariantCaller::IsGoodAltAlleleWithReason(
    const Allele& allele, const int total_count,
    const bool apply_trio_coefficient) const {
  if (allele.type() == AlleleType::REFERENCE) {
    return AlleleRejectionAcceptance::REJECTED_REF;
  }

  if (allele.count() < min_count(allele)) {
    return AlleleRejectionAcceptance::REJECTED_LOW_SUPPORT;
  }

  if (allele.type() == AlleleType::SOFT_CLIP) {
    return AlleleRejectionAcceptance::REJECTED_OTHER;
  }

  if ((1.0 * allele.count()) / total_count <
      min_fraction(allele) *
          (apply_trio_coefficient ? options_.min_fraction_multiplier() : 1.0)) {
    return AlleleRejectionAcceptance::REJECTED_LOW_RATIO;
  }

  return VariantCaller::AlleleRejectionAcceptance::ACCEPTED;
}

bool IsAllelesTheSame(const Allele& allele1, const Allele& allele2) {
  return (allele1.bases() == allele2.bases() &&
          allele1.type() == allele2.type());
}

// Select the subset of GoodAltAlleles from the alleles of allele_count.
//
// Returns the vector of allele objects from allele_count that satisfy
// IsGoodAltAllele().
std::vector<Allele> VariantCaller::SelectAltAlleles(
    const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
    absl::string_view target_sample) const {
  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  const AlleleCount& target_sample_allele_count =
      allele_counts.at(target_sample);
  std::vector<AlleleCount> all_samples_allele_counts;
  all_samples_allele_counts.reserve(allele_counts.size());
  // "Non-target" samples are referring to all the samples that are providing
  // supportive information. Usually the main truth labels are not from this
  // sample, or usually it means that the calls coming from these non-target
  // samples are not the main focus of our problem.
  std::vector<AlleleCount> non_target_allele_counts;
  non_target_allele_counts.reserve(allele_counts.size());
  for (const auto& allele_counts_entry : allele_counts) {
    all_samples_allele_counts.push_back(allele_counts_entry.second);
    if (allele_counts_entry.first != target_sample) {
      non_target_allele_counts.push_back(allele_counts_entry.second);
    }
  }

  const std::vector<Allele> target_sample_alleles =
      SumAlleleCounts(target_sample_allele_count);
  const std::vector<Allele> all_sample_alleles =
      SumAlleleCounts(all_samples_allele_counts);
  const std::vector<Allele> non_target_sample_alleles =
      SumAlleleCounts(non_target_allele_counts);

  const int target_samples_total_count =
      TotalAlleleCounts(target_sample_allele_count);
  const int all_samples_total_count =
      TotalAlleleCounts(all_samples_allele_counts);

  std::vector<Allele> alt_alleles;
  // First process target_sample_alleles
  for (const auto& allele : target_sample_alleles) {
    bool skip_high_af_allele_for_non_target = false;
    // Having a double for-loop seems inefficient. Can be room for improvement.
    for (const auto& non_target_sample_allele : non_target_sample_alleles) {
      if (!IsAllelesTheSame(allele, non_target_sample_allele)) continue;
      int non_target_total_count =
          all_samples_total_count - target_samples_total_count;
      float max_fraction_for_non_target_sample =
          non_target_sample_allele.type() == AlleleType::SUBSTITUTION
          ? options_.max_fraction_snps_for_non_target_sample()
          : options_.max_fraction_indels_for_non_target_sample();
      if (max_fraction_for_non_target_sample > 0 &&
          (1.0 * non_target_sample_allele.count() / non_target_total_count) >
          max_fraction_for_non_target_sample) {
        skip_high_af_allele_for_non_target = true;
        break;
      }
    }
    if (skip_high_af_allele_for_non_target) continue;

    AlleleRejectionAcceptance allele_acceptance =
        IsGoodAltAlleleWithReason(allele, target_samples_total_count, false);
    if (allele_acceptance == AlleleRejectionAcceptance::ACCEPTED) {
      alt_alleles.push_back(allele);
      continue;
    }
    if (allele_acceptance == AlleleRejectionAcceptance::REJECTED_LOW_RATIO ||
        allele_acceptance == AlleleRejectionAcceptance::REJECTED_LOW_SUPPORT) {
      for (const auto& all_samples_allele : all_sample_alleles) {
        if (IsAllelesTheSame(allele, all_samples_allele) &&
            AlleleRejectionAcceptance::ACCEPTED ==
                IsGoodAltAlleleWithReason(all_samples_allele,
                                          all_samples_total_count,
                                true)) {
          alt_alleles.push_back(allele);
          break;
        }  // if (Found good allele in other samples)
      }    // for (all alleles in all samples)
    }      // if (allele rejected)
  }        // for (alleles in target samples)

  return alt_alleles;
}

// Adds a single VariantCall with sample_name, genotypes, and gq (bound to the
// "GQ" key of info with a numerical value of gq, if provided) to variant.
void AddGenotypes(const std::string& sample_name,
                  const std::vector<int>& genotypes, Variant* variant) {
  CHECK(variant != nullptr);

  VariantCall* call = variant->add_calls();
  call->set_call_set_name(sample_name);
  for (const auto genotype : genotypes) {
    call->add_genotype(genotype);
  }
}

AlleleMap BuildAlleleMap(const AlleleCount& allele_count,
                         const std::vector<Allele>& alt_alleles,
                         const std::string& ref_bases) {
  AlleleMap allele_map;

  // Compute the alt alleles, recording the mapping from each Allele to its
  // corresponding allele in the Variant format.
  for (const auto& alt_allele : alt_alleles) {
    const std::string_view alt_bases = alt_allele.bases();
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

AlleleMap RemoveInvalidDels(const AlleleMap& allele_map,
                            absl::string_view ref_bases) {
  AlleleMap allele_map_mod;
  absl::btree_map<Allele, int, OrderAllele> read_counts;
  int num_of_dels = 0;
  bool has_deletion_adjacent_to_snp = false;

  // Search for deletions and check if there is a deletion with the preceding
  // SNP. SNP is followed by deletion if deletion's alt base is different from
  // the ref.
  // In addition, read count is stored for each deletion.
  for (const auto& elt : allele_map) {
    if (elt.first.type() == AlleleType::DELETION) {
      read_counts[elt.first] += elt.first.count();
      num_of_dels++;
      if (elt.second[0] != ref_bases[0]) {
        has_deletion_adjacent_to_snp = true;
      }
    }
  }

  // If more than 1 DELs and their alt bases are different we need to keep just
  // one. The one with higher read support is kept.
  if (num_of_dels > 1 && has_deletion_adjacent_to_snp) {
    Allele max_allele =
        std::max_element(read_counts.begin(), read_counts.end(),
                         [](const std::pair<const Allele, int>& element1,
                            const std::pair<const Allele, int>& element2) {
                           return element1.second < element2.second;
                         })
            ->first;

    if (!max_allele.bases().empty()) {
      for (const auto& elt : allele_map) {
        if ((elt.first.type() == AlleleType::DELETION &&
             elt.second == allele_map.at(max_allele)) ||
            elt.first.type() != AlleleType::DELETION) {
          allele_map_mod[elt.first] = elt.second;
        }
      }
      return allele_map_mod;
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

    absl::btree_map<absl::string_view, const Allele*> alt_to_alleles;
    for (const auto& entry : allele_map) {
      alt_to_alleles[entry.second] = &entry.first;
    }
    CHECK(alt_to_alleles.size() == allele_map.size())
        << "Non-unique alternative alleles!";
    for (const std::string& alt : variant->alternate_bases()) {
      const Allele& allele = *alt_to_alleles.find(alt)->second;
      ad.push_back(allele.count());
      vaf.push_back(1.0 * allele.count() / dp);
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

template <class T>
std::vector<T> VariantCaller::AlleleCountsGenerator(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
    const std::string& target_sample,
    std::optional<T> (VariantCaller::*F)(
        const absl::node_hash_map<std::string, AlleleCount>&,
        const std::string&, const std::vector<AlleleCount>*,
        std::vector<AlleleCount>::const_iterator*) const) const {
  // Get Allele counts for the target sample
  auto it = allele_counters.find(target_sample);
  if (it == allele_counters.end()) {
    LOG(WARNING)
        << "allele_counters collection does not contain target sample!";
    return std::vector<T>();
  }

  // Contains AlleleCount objects for each position of the target sample.
  const std::vector<AlleleCount>& target_sample_allele_counts =
      it->second->Counts();

  // Initialize a vector of iterators - one iterator per sample.
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;
  for (const auto& sample_allele_counters : allele_counters) {
    allele_counter_iterators[sample_allele_counters.first] =
        sample_allele_counters.second->Counts().begin();
  }

  std::vector<T> items;

  // Iterate through AlleleCount objects for each position, moving iterators
  // for each sample simultaneously.
  while (allele_counter_iterators[target_sample] !=
         target_sample_allele_counts.end()) {
    absl::node_hash_map<std::string, AlleleCount> allele_counts_per_sample;
    for (const auto& sample_counter : allele_counters) {
      if (allele_counter_iterators[sample_counter.first] !=
          allele_counters.at(sample_counter.first)->Counts().end()) {
        // allele_counts_per_sample contain AlleleCount for each sample for one
        // position.
        allele_counts_per_sample[sample_counter.first] =
            *(allele_counter_iterators[sample_counter.first]);
      }
    }
    // Calling CallVariant for one position. allele_counts_per_sample contains
    // AlleleCount object for this position for each sample.
    std::optional<T> item = (this->*F)(
        allele_counts_per_sample, target_sample, &target_sample_allele_counts,
        &allele_counter_iterators[target_sample]);
    if (item) {
      items.push_back(*item);
    }

    // Increment all iterators.
    for (auto& it : allele_counter_iterators) {
      if (it.second != allele_counters.at(it.first)->Counts().end()) {
        it.second++;
      }
    }
  }
  return items;
}

std::vector<DeepVariantCall> VariantCaller::CallsFromAlleleCounts(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
    const std::string& target_sample) const {
  return AlleleCountsGenerator<DeepVariantCall>(allele_counters, target_sample,
                                                &VariantCaller::CallVariant);
}

std::vector<int> VariantCaller::CallPositionsFromAlleleCounts(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
    const std::string& target_sample) const {
  return AlleleCountsGenerator<int>(allele_counters, target_sample,
                                    &VariantCaller::CallVariantPosition);
}

std::optional<int> VariantCaller::CallVariantPosition(
    const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
    const std::string& target_sample,
    const std::vector<AlleleCount>* target_sample_allele_counts,
    std::vector<AlleleCount>::const_iterator*
        target_sample_allele_count_iterator) const {
  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  const AlleleCount& target_sample_allele_count =
      allele_counts.at(target_sample);
  if (!nucleus::AreCanonicalBases(target_sample_allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return std::nullopt;
  }

  const std::vector<Allele> alt_alleles =
      SelectAltAlleles(allele_counts, target_sample);
  if (alt_alleles.empty() && !KeepReferenceSite()) {
    return std::nullopt;
  }
  return target_sample_allele_count.position().position();
}

std::optional<DeepVariantCall> VariantCaller::CallVariant(
    const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
    const std::string& target_sample,
    const std::vector<AlleleCount>* target_sample_allele_counts,
    std::vector<AlleleCount>::const_iterator*
        target_sample_allele_count_iterator) const {
  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  const AlleleCount& target_sample_allele_count =
      allele_counts.at(target_sample);
  if (!nucleus::AreCanonicalBases(target_sample_allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return std::nullopt;
  }

  const std::vector<Allele> alt_alleles =
      SelectAltAlleles(allele_counts, target_sample);
  if (alt_alleles.empty() && !KeepReferenceSite()) {
    return std::nullopt;
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
  variant->set_reference_name(
      target_sample_allele_count.position().reference_name());
  variant->set_start(target_sample_allele_count.position().position());
  const std::string refbases =
      CalcRefBases(target_sample_allele_count.ref_base(), alt_alleles);
  variant->set_reference_bases(refbases);
  variant->set_end(variant->start() + refbases.size());
  AddGenotypes(options_.sample_name(), {-1, -1}, variant);

  // Compute the map from read alleles to the alleles we'll use in our Variant.
  // Add the alternate alleles from our allele_map to the variant.
  const AlleleMap allele_map =
      BuildAlleleMap(target_sample_allele_count, alt_alleles, refbases);
  for (const auto& elt : allele_map) {
    variant->add_alternate_bases(elt.second);
  }
  // If we don't have any alt_alleles, we are generating a reference site so
  // add in the kNoAltAllele.
  if (alt_alleles.empty()) variant->add_alternate_bases(kNoAltAllele);
  std::sort(variant->mutable_alternate_bases()->pointer_begin(),
            variant->mutable_alternate_bases()->pointer_end(),
            StringPtrLessThan());

  AddReadDepths(target_sample_allele_count, allele_map, variant);
  AddSupportingReads(allele_counts, allele_map, target_sample, &call);
  if (options_.small_model_vaf_context_window_size() > 0) {
    AddAdjacentAlleleFractionsAtPosition(
        options_.small_model_vaf_context_window_size(),
        *target_sample_allele_counts, *target_sample_allele_count_iterator,
        &call);
  }

  return std::make_optional(call);
}

AlleleMap::const_iterator FindAllele(const Allele& allele,
                                     const AlleleMap& allele_map) {
  for (auto it = allele_map.begin(); it != allele_map.end(); it++) {
    if (IsAllelesTheSame(it->first, allele)) {
      return it;
    }
  }
  return allele_map.end();
}

void VariantCaller::AddSupportingReads(
    const absl::node_hash_map<std::string, AlleleCount>& allele_counts,
    const AlleleMap& allele_map, const std::string& target_sample,
    DeepVariantCall* call) const {
  // Iterate over each read in the allele_count, and add its name to the
  // supporting reads of for the Variant allele it supports.
  const std::string unknown_allele = kSupportingUncalledAllele;
  absl::flat_hash_map<std::string, absl::flat_hash_set<std::string>>
      alt_allele_support;
  absl::flat_hash_set<std::string> ref_support;
  for (const auto& allele_counts_entry : allele_counts) {
    const AlleleCount& allele_count = allele_counts_entry.second;
    for (const auto& read_name_allele : allele_count.read_alleles()) {
      const std::string& read_name = read_name_allele.first;
      const Allele& allele = read_name_allele.second;

      // Skip reference supporting reads, as they aren't included in the
      // supporting reads for alternate alleles.
      if (allele.type() != AlleleType::REFERENCE) {
        auto it = FindAllele(allele, allele_map);
        const std::string& supported_allele =
            it == allele_map.end() ? unknown_allele : it->second;
        DeepVariantCall::SupportingReads& supports =
            (*call->mutable_allele_support())[supported_allele];
        // Check that this read does not exist in supports already. It may
        // happen if candidate is created from multiple samples and read with
        // the same id exists in multiple samples. Multiple problems may arise
        // from this: number of supporting reads is calculated incorrectly,
        // phasing may not work due to loops in the graph caused by multiple
        // reads with the same id supporting the same allele.
        auto [new_item, is_inserted] =
            alt_allele_support[supported_allele].insert(read_name);
        if (!is_inserted) continue;

        supports.add_read_names(read_name);
        DeepVariantCall_SupportingReadsExt& support_infos =
            (*call->mutable_allele_support_ext())[supported_allele];
        DeepVariantCall_ReadSupport* read_info = support_infos.add_read_infos();
        read_info->set_read_name(read_name);
        read_info->set_is_low_quality(allele.is_low_quality());
        read_info->set_mapping_quality(allele.mapping_quality());
        read_info->set_average_base_quality(allele.avg_base_quality());
      } else {
        call->add_ref_support(read_name);
        DeepVariantCall_SupportingReadsExt& support_infos =
            (*call->mutable_ref_support_ext());
        DeepVariantCall_ReadSupport* read_info = support_infos.add_read_infos();
        auto [new_element, is_inserted] = ref_support.insert(read_name);
        if (!is_inserted) continue;

        read_info->set_read_name(read_name);
        read_info->set_is_low_quality(allele.is_low_quality());
        read_info->set_mapping_quality(allele.mapping_quality());
        read_info->set_average_base_quality(allele.avg_base_quality());
      }
    }
  }
}

void VariantCaller::AddAdjacentAlleleFractionsAtPosition(
    const int window_size,
    const std::vector<AlleleCount>& target_sample_allele_counts,
    std::vector<AlleleCount>::const_iterator
        target_sample_allele_count_iterator,
    DeepVariantCall* call) const {
  int index = std::distance(target_sample_allele_counts.begin(),
                            target_sample_allele_count_iterator);
  int half_window_size = window_size / 2;
  int start = std::min(index, half_window_size);
  int end = std::min((int)target_sample_allele_counts.size() - index,
                     half_window_size + 1);
  const std::vector<AlleleCount>& context_allele_counts = {
      target_sample_allele_count_iterator - start,
      target_sample_allele_count_iterator + end};
  for (const auto& context_allele_count : context_allele_counts) {
    int depth = context_allele_count.ref_supporting_read_count() +
                context_allele_count.read_alleles_size();
    int vaf = 0;
    if (depth > 0) {
      vaf = (100 * context_allele_count.read_alleles_size()) / depth;
    }
    (*call->mutable_allele_frequency_at_position())
        [context_allele_count.position().position()] = vaf;
  }
}

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
