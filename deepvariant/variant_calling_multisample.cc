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
#include <cstdint>
#include <iterator>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/variants.pb.h"
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
const char* const kMFFormatField = "MF";
const char* const kMDFormatField = "MD";

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
                         absl::Span<const Allele> alt_alleles) {
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
                          absl::string_view variant_ref, const uint32_t from) {
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

bool VariantCaller::AlleleFilter(
    const Allele& allele, const AlleleCount& target_sample_allele_count,
    absl::Span<const AlleleCount> all_samples_allele_counts,
    absl::Span<const Allele> non_target_sample_alleles) const {
  const int target_samples_total_count =
      TotalAlleleCounts(target_sample_allele_count);
  const int all_samples_total_count =
      TotalAlleleCounts(all_samples_allele_counts);
  const std::vector<Allele> all_sample_alleles =
      SumAlleleCounts(all_samples_allele_counts);

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
  if (skip_high_af_allele_for_non_target) return false;

  AlleleRejectionAcceptance allele_acceptance =
      IsGoodAltAlleleWithReason(allele, target_samples_total_count, false);
  if (allele_acceptance == AlleleRejectionAcceptance::ACCEPTED) {
    return true;
  }
  if (allele_acceptance == AlleleRejectionAcceptance::REJECTED_LOW_RATIO ||
      allele_acceptance == AlleleRejectionAcceptance::REJECTED_LOW_SUPPORT) {
    for (const auto& all_samples_allele : all_sample_alleles) {
      if (IsAllelesTheSame(allele, all_samples_allele) &&
          AlleleRejectionAcceptance::ACCEPTED ==
              IsGoodAltAlleleWithReason(all_samples_allele,
                                        all_samples_total_count,
                              true)) {
        return true;
      }  // if (Found good allele in other samples)
    }    // for (all alleles in all samples)
  }      // if (allele rejected)
  return false;
}

namespace {
  // Helper function to generate a map from read id to a vector of alt alleles
  // that are supported by each read.
  absl::flat_hash_map<
      std::string, std::vector<AlleleAtPosition>> CreateCombinedAllelesSupport(
        const std::vector<AlleleCount>& allele_counts_context,
        absl::string_view del_allele_ref_bases,
        int del_start,
        int del_len
      ) {
    absl::flat_hash_map<std::string, std::vector<AlleleAtPosition>>
        read_to_alt_alleles;
    int found_alt_allele_overlapped_by_deletion = 0;
    bool overlapping_del_found = false;
    // Iterating over all positions starting at the del_start.
    for (const auto& allele_count : allele_counts_context) {
      int allele_pos = allele_count.position().position();
      if (allele_pos < del_start) {
        continue;
      }
      if (allele_pos >= del_start + del_len) {
        break;
      }
      for (const auto& [read_id, read_allele] :
            allele_count.read_alleles()) {
        // Skip alleles for the deletion itself.
        if (allele_pos == del_start &&
            read_allele.type() == AlleleType::DELETION &&
            read_allele.bases().size() == del_len) {
          continue;
        }
        // We cannot create complex variant if there are other deletions
        // overlapping our deletion.
        if (read_allele.type() == AlleleType::DELETION) {
          found_alt_allele_overlapped_by_deletion = 0;
          overlapping_del_found = true;
          break;
        }
        // We are only interested in cases where there is alt allele that
        // starts after the deletion start and ends before the deletion
        // end.
        if (allele_pos >= del_start &&
            read_allele.type() != AlleleType::REFERENCE) {
          found_alt_allele_overlapped_by_deletion++;
        }
        read_to_alt_alleles[read_id].push_back(
            {.alt_bases = read_allele.bases(),
              .type = read_allele.type(),
              .position = allele_pos});
      }  // for (read_id, read_allele)
    }  // for (allele_counts_context)
    if (found_alt_allele_overlapped_by_deletion < 1 || overlapping_del_found) {
      read_to_alt_alleles.clear();
    }
    return read_to_alt_alleles;
  }
}  // namespace


// Helper function to generate a map of complex alleles to supporting reads.
// This function creates complex alleles by concatenating alt alleles that are
// supported by the same read. Results are returned as a map from complex
// allele to a vector of supporting reads. It is expected that all alt alleles
// are not deletions.
// Args:
//   read_to_alt_alleles: map from read id to a vector of alt alleles that are
//     supported by the read.
//   del_start: start position of the deletion.
//   del_len: length of the deletion.
//   del_allele_ref_bases: reference bases of the deletion.
absl::flat_hash_map<std::string, std::vector<std::string>>
CreateComplexAllelesSupport(
    const absl::flat_hash_map<std::string, std::vector<AlleleAtPosition>>&
        read_to_alt_alleles,
    int del_start, int del_len, absl::string_view del_allele_ref_bases) {
  absl::flat_hash_map<std::string, std::vector<std::string>>
      complex_allele_to_reads;
  // Iterating over all potential complex alleles supported by a read.
  for (const auto& [read_id, alt_alleles] : read_to_alt_alleles) {
    int start_pos = 0;
    // The set is used to check the uniqueness of complex alleles.
    absl::flat_hash_set<
        std::tuple<AlleleType, std::string>> complex_alleles_strings;
    std::string complex_allele;
    // Iterate over alt alleles for this read.
    for (const auto& allele : alt_alleles) {
      if (allele.type == AlleleType::UNSPECIFIED) {
        continue;
      }
      // Deletion alleles should not be passed to this function.
      CHECK_NE(allele.type, AlleleType::DELETION);
      const int relative_alt_allele_pos = allele.position - del_start;
      // Add leading ref bases.
      if (relative_alt_allele_pos > start_pos &&
          relative_alt_allele_pos <= del_len) {
        absl::StrAppend(&complex_allele,
                        del_allele_ref_bases.substr(
                            start_pos, relative_alt_allele_pos - start_pos));
        start_pos += relative_alt_allele_pos - start_pos;
      }
      // Add allele bases.
      absl::StrAppend(&complex_allele, allele.alt_bases);
      // Update the start_pos.
      if (allele.type != AlleleType::INSERTION) {
        start_pos = relative_alt_allele_pos + allele.alt_bases.size();
      } else {
        start_pos += 1;
      }
    }
    // Record the complex allele if it is not empty and it is not a duplicate.
    auto complex_allele_and_type = std::make_tuple(
        AlleleTypeFromAlt(del_allele_ref_bases, complex_allele),
        complex_allele);
    if (!complex_allele.empty() && start_pos <= del_len
        && !complex_alleles_strings.contains(complex_allele_and_type)) {
      // Add trailing ref bases.
      absl::StrAppend(&complex_allele,
                      del_allele_ref_bases.substr(start_pos));
      complex_allele_to_reads[complex_allele].push_back(read_id);
      complex_alleles_strings.insert(complex_allele_and_type);
    } else {
      // For now we drop sites where at least one of the complex alleles
      // couldn't be generated. In this case normal alleles will be used.
      return {};
    }
  }
  return complex_allele_to_reads;
}


void ReassignReadSupportForComplexAlleles(
    const absl::node_hash_map<std::string, AlleleCount>& original_allele_counts,
    absl::string_view target_sample,
    const absl::flat_hash_map<std::string, std::vector<std::string>>&
      complex_allele_to_reads,
    absl::string_view del_allele_ref_bases,
    absl::node_hash_map<std::string, AlleleCount>& allele_counts_mod) {
  bool alt_allels_were_modified = false;
  AlleleCount target_sample_allele_count_mod;
  auto it  = original_allele_counts.find(target_sample);
  CHECK(it != original_allele_counts.end());
  const AlleleCount& target_sample_allele_count = it->second;
  target_sample_allele_count_mod = target_sample_allele_count;
  int total_ref_count = 0;
  for (const auto& [complex_allele, reads] : complex_allele_to_reads) {
    alt_allels_were_modified = true;
    for (const auto& read_id : reads) {
      // TODO: With this check we loose supporting reads that start
      // after the deletion's starting position. This implementation doesn't
      // handle such cases. In the future we may want to add the functionality
      // that would allow to have two reads to support a complex allele.
      // Complex allele: ACGTCTATG
      //                 |||||||||
      // Read1:      ACTGACGT|||||
      // Read2:              CTATGATC
      // It is not a problem for long reads since it would be a rare case when
      // reads break in the middle of the complex allele. But, for short reads
      // it is a problem.
      if (!target_sample_allele_count_mod.read_alleles().contains(read_id))  {
        continue;
      }
      auto& read_allele =
          target_sample_allele_count_mod.mutable_read_alleles()->at(read_id);
      if (complex_allele == del_allele_ref_bases) {
        read_allele.set_type(AlleleType::REFERENCE);
        total_ref_count++;
      } else {
        read_allele.set_type(AlleleType::SUBSTITUTION);
      }
      read_allele.set_bases(complex_allele);
    }
  }
  target_sample_allele_count_mod.set_ref_supporting_read_count(total_ref_count);
  // Update allele counts for all samples. Currently we only update the target
  // sample.
  allele_counts_mod.clear();
  for (const auto& [sample, allele_count] : original_allele_counts) {
    if (sample == target_sample) {
      allele_counts_mod.insert({sample, target_sample_allele_count_mod});
    } else {
      allele_counts_mod.insert({sample, allele_count});
    }
  }
}

// Implementation of complex variant representation.
//
// This function is called when we have a deletion allele in the target sample.
// We need to check if there are other alleles that overlap with this deletion.
// If there are, we need to create a complex variant allele that is a
// concatenation of all the overlapping alleles. In this case complex allele is
// an allele that is created by concatenating all simple alleles that are
// overlapped by the deletion. Complex allele is a combination of SNP, REF,
// and INS.
// For example, if we have a deletion at position 10 and SNP 1 at position 10,
// and SNP 2 at position 12. And SNP 1 and SNP 2 are supported by the same set
// of reads.
//   Position 10 Allele 1: ATCG -> A (Deletion)
//   Position 10 Allele 2: A -> T    (SNP 1)
//   Position 12 Allele 3: C -> A    (SNP 2)
//
// This function will create a complex variant:
//   Position 10 Allele 1: ATCG -> A       (Deletion)
//   Position 10 Allele 2: ATCG -> TTAG    (SNP 1, 2)
SelectAltAllelesResult
  VariantCaller::SelectAltAllelesWithComplexVariant(
        const SelectAltAllelesWithComplexVariantInputOptions& options) const {
  SelectAltAllelesResult output_options;
  // Scan alt_alleles for deletion and create complex variant candidate if there
  // are other alleles overlapped with this deletion.
  std::vector<Allele>::const_iterator allele_with_del = std::find_if(
      options.alt_alleles.cbegin(), options.alt_alleles.cend(),
      [](const Allele& allele) {
        return allele.type() == AlleleType::DELETION; });
  if (allele_with_del == options.alt_alleles.cend()) {
    return {
      .alt_alleles = options.alt_alleles,
      .complex_variant_created = false,
      .ref_bases = "",
    };
  }

  auto it = allele_counters_per_sample_.find(target_sample_);
  CHECK(it != allele_counters_per_sample_.end());
  const std::vector<AlleleCount>& allele_counts_context = it->second->Counts();

  auto it_allele_count = options.allele_counts_by_sample.find(target_sample_);
  CHECK(it_allele_count != options.allele_counts_by_sample.end());
  const AlleleCount& target_sample_allele_count = it_allele_count->second;
  // Since there may be multiple deletions of different lengths we need to
  // set del_len to the largest one. This is done in CalcRefBases().
  output_options.ref_bases =
      CalcRefBases(target_sample_allele_count.ref_base(), options.alt_alleles);
  // If there are multiple deletions we take the largest one.
  int del_len = output_options.ref_bases.size();
  int del_start = target_sample_allele_count.position().position();

  absl::flat_hash_map<std::string, std::vector<AlleleAtPosition>>
      read_to_alt_alleles = CreateCombinedAllelesSupport(
          allele_counts_context, output_options.ref_bases,
          del_start,  // Starting position of the deletion allele.
          del_len);  // Length of the deletion allele. This will determine the
                     // length of the complex allele.
  if (read_to_alt_alleles.empty()) {
    output_options.complex_variant_created = false;
    output_options.ref_bases = "";
  } else {
    output_options.complex_variant_created = true;
  }

  absl::flat_hash_map<std::string, std::vector<std::string>>
      complex_allele_to_reads = CreateComplexAllelesSupport(
          read_to_alt_alleles, del_start, del_len, output_options.ref_bases);

  ReassignReadSupportForComplexAlleles(
      options.allele_counts_by_sample,
      target_sample_,
      complex_allele_to_reads,
      output_options.ref_bases,
      output_options.allele_counts_mod);

  // Call SelectAltAlleles one more time with the modified allele counts to run
  // newly created complex alleles through the filtering.
  SelectAltAllelesResult output_options_temp = SelectAltAlleles(
    {
      .allele_counts_by_sample = output_options.allele_counts_mod,
      .create_complex_alleles = false,
      .prev_deletion_end = 0,
    }
  );
  output_options.alt_alleles = std::move(output_options_temp).alt_alleles;
  return output_options;
}

// Select the subset of GoodAltAlleles from the alleles of allele_count at a
// single position.
//
// Returns the vector of allele objects from allele_count that satisfy
// IsGoodAltAllele().
SelectAltAllelesResult VariantCaller::SelectAltAlleles(
    const SelectAltAllelesInputOptions& options) const {
  SelectAltAllelesResult output_options;
  // allele_counts_mod is initialized to be the same as allele_counts. If
  // complex variants are eneabled then allele_counts_mod will be modified.
  for (const auto& [sample_name, allele_count] :
           options.allele_counts_by_sample) {
      output_options.allele_counts_mod.insert({sample_name, allele_count});
  }

  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  auto it = options.allele_counts_by_sample.find(target_sample_);
  CHECK(it != options.allele_counts_by_sample.end());
  const AlleleCount& target_sample_allele_count = it->second;
  std::vector<AlleleCount> all_samples_allele_counts;
  all_samples_allele_counts.reserve(options.allele_counts_by_sample.size());
  // "Non-target" samples are referring to all the samples that are providing
  // supportive information. Usually the main truth labels are not from this
  // sample, or usually it means that the calls coming from these non-target
  // samples are not the main focus of our problem.
  std::vector<AlleleCount> non_target_allele_counts;
  non_target_allele_counts.reserve(options.allele_counts_by_sample.size());
  for (const auto& allele_counts_entry : options.allele_counts_by_sample) {
    all_samples_allele_counts.push_back(allele_counts_entry.second);
    if (allele_counts_entry.first != target_sample_) {
      non_target_allele_counts.push_back(allele_counts_entry.second);
    }
  }

  const std::vector<Allele> target_sample_alleles =
      SumAlleleCounts(target_sample_allele_count);
  const std::vector<Allele> all_sample_alleles =
      SumAlleleCounts(all_samples_allele_counts);
  const std::vector<Allele> non_target_sample_alleles =
      SumAlleleCounts(non_target_allele_counts);

  std::vector<Allele> alt_alleles;
  for (const auto& allele : target_sample_alleles) {
    if (AlleleFilter(
            allele, target_sample_allele_count, all_samples_allele_counts,
            non_target_sample_alleles)) {
      alt_alleles.push_back(allele);
    }
  }  // for (alleles in target samples)

  if (options.create_complex_alleles &&
      options.prev_deletion_end <=
      target_sample_allele_count.position().position()) {
    return SelectAltAllelesWithComplexVariant(
        {
          .allele_counts_by_sample = options.allele_counts_by_sample,
          .alt_alleles = std::move(alt_alleles)
        });
  } else {
    return {
      .alt_alleles = alt_alleles,
      .complex_variant_created = false,
      .ref_bases = "",
    };
  }
}

// Adds a single VariantCall with sample_name, genotypes, and gq (bound to the
// "GQ" key of info with a numerical value of gq, if provided) to variant.
void AddGenotypes(const std::string& sample_name,
                  absl::Span<const int> genotypes, Variant* variant) {
  CHECK(variant != nullptr);

  VariantCall* call = variant->add_calls();
  call->set_call_set_name(sample_name);
  for (const auto genotype : genotypes) {
    call->add_genotype(genotype);
  }
}

AlleleMap BuildAlleleMap(const AlleleCount& allele_count,
                         absl::Span<const Allele> alt_alleles,
                         absl::string_view ref_bases) {
  AlleleMap allele_map;

  // Compute the alt alleles, recording the mapping from each Allele to its
  // corresponding allele in the Variant format.
  for (const auto& alt_allele : alt_alleles) {
    const std::string_view alt_bases = alt_allele.bases();
    switch (alt_allele.type()) {
      case AlleleType::SUBSTITUTION:
        if (alt_bases.size() > 1 && ref_bases.size() > 1) {
          allele_map[alt_allele] = alt_bases;
        } else {
          allele_map[alt_allele] = MakeAltAllele(alt_bases, ref_bases, 1);
        }
        break;
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
    for (const auto& [allele, alt_bases] : allele_map) {
      alt_to_alleles[alt_bases] = &allele;
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

// Returns true if the current site should be emitted, even if it's a reference
// site. This function is used to return reference site samples if the
// member variable fraction_reference_sites_to_emit >= 0.0 by pulling draws
// from a random number and returning true if the number <= the threshold.
bool VariantCaller::KeepReferenceSite() const {
  return options_.fraction_reference_sites_to_emit() > 0.0 && sampler_.Keep();
}

template <class T>
std::vector<T> VariantCaller::AlleleCountsGenerator(
    std::optional<T> (VariantCaller::*F)(
        const absl::node_hash_map<std::string,
        AlleleCount>&,
        std::vector<AlleleCount>::const_iterator*,
        int& skip_next_count_param,
        int& prev_deletion_end_param) const) const {
  // Initialize a vector of iterators - one iterator per sample.
  absl::node_hash_map<std::string, std::vector<AlleleCount>::const_iterator>
      allele_counter_iterators;
  for (const auto& [sample, allele_counter] : allele_counters_per_sample_) {
    allele_counter_iterators[sample] =
        allele_counter->Counts().begin();
  }
  const std::vector<AlleleCount>& target_sample_allele_counts =
      allele_counters_per_sample_.at(target_sample_)->Counts();

  std::vector<T> items;

  // Iterate through AlleleCount objects for each position, moving iterators
  // for each sample simultaneously.
  int skip_next_count = 0;
  int prev_deletion_end = 0;
  while (allele_counter_iterators[target_sample_] !=
         target_sample_allele_counts.end()) {
    absl::node_hash_map<std::string, AlleleCount> allele_counts_per_sample;
    for (const auto& sample_counter : allele_counters_per_sample_) {
      if (allele_counter_iterators[sample_counter.first] !=
          allele_counters_per_sample_.at(
              sample_counter.first)->Counts().end()) {
        // allele_counts_per_sample contain AlleleCount for each sample for one
        // position.
        allele_counts_per_sample[sample_counter.first] =
            *(allele_counter_iterators[sample_counter.first]);
      }
    }
    // Calling CallVariant for one position. allele_counts_per_sample contains
    // AlleleCount object for this position for each sample.
    if (skip_next_count > 0) {
      skip_next_count--;
    } else {
      std::optional<T> item = (this->*F)(
          allele_counts_per_sample,
          &allele_counter_iterators[target_sample_],
          skip_next_count,
          prev_deletion_end);
      if (item) {
        items.push_back(*item);
      }
    }

    // Increment all iterators.
    for (auto& it : allele_counter_iterators) {
      if (it.second !=
          allele_counters_per_sample_.at(it.first)->Counts().end()) {
        it.second++;
      }
    }
  }
  return items;
}

std::vector<DeepVariantCall> VariantCaller::CallsFromAlleleCounts(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
    const std::string& target_sample) {
  this->allele_counters_per_sample_ = allele_counters;
  this->target_sample_ = target_sample;
  // Get Allele counts for the target sample
  auto it = allele_counters.find(target_sample);
  if (it == allele_counters.end()) {
    LOG(FATAL)
        << "allele_counters collection does not contain target sample!";
  }
  // Merge methylated supporting reads in the negative strand
  // (Pacbio 5mC specific)
  if (options_.enable_methylation_aware_phasing()) {
    MergeMethylatedAlleleCounts(allele_counters_per_sample_);
  }
  return AlleleCountsGenerator<DeepVariantCall>(&VariantCaller::CallVariant);
}

std::vector<int> VariantCaller::CallPositionsFromAlleleCounts(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters,
    const std::string& target_sample) {
  this->allele_counters_per_sample_ = allele_counters;
  this->target_sample_ = target_sample;
  // Get Allele counts for the target sample
  auto it = allele_counters.find(target_sample);
  if (it == allele_counters.end()) {
    LOG(FATAL)
        << "allele_counters collection does not contain target sample!";
  }
  return AlleleCountsGenerator<int>(&VariantCaller::CallVariantPosition);
}

std::optional<int> VariantCaller::CallVariantPosition(
    const absl::node_hash_map<std::string, AlleleCount>&
        allele_counts_by_sample,
    std::vector<AlleleCount>::const_iterator*
        target_sample_allele_count_iterator,
    int& skip_next_count,
    int& prev_deletion_end) const {
  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  const AlleleCount& target_sample_allele_count =
      allele_counts_by_sample.at(target_sample_);
  if (!nucleus::AreCanonicalBases(target_sample_allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return std::nullopt;
  }
  SelectAltAllelesResult output_options = SelectAltAlleles(
      {
      .allele_counts_by_sample = allele_counts_by_sample,
      .create_complex_alleles = false,
      .prev_deletion_end = prev_deletion_end
    });

  // Include reference site as candidate for methylation-aware phasing if it is
  // methylated.
  // However, if the methylated reference site is in X or Y chromosome,
  // we do not include it as a candidate.
  bool has_methylation = false;
  std::string chrom = target_sample_allele_count.position().reference_name();
  if (options_.enable_methylation_aware_phasing() &&
      IsReferenceSite(target_sample_allele_count) &&
      !IsExcludedMethylationContig(chrom)) {
    // Count methylated reads in reference sites
    for (const auto& read_entry : target_sample_allele_count.read_alleles()) {
      const Allele& allele = read_entry.second;
      if (allele.is_methylated()) {
        has_methylation = true;
        break;
      }
    }
  }

  if (!KeepReferenceSite() && output_options.alt_alleles.empty() &&
      !has_methylation) {
    return std::nullopt;
  }

  return target_sample_allele_count.position().position();
}

std::optional<DeepVariantCall> VariantCaller::CallVariant(
    const absl::node_hash_map<std::string, AlleleCount>&
        allele_counts_per_sample,
    std::vector<AlleleCount>::const_iterator*
        target_sample_allele_count_iterator,
    int& skip_next_count,
    int& prev_deletion_end) const {
  // allele_counts.at will throw an exception if key is not found.
  // Absent target_sample is a critical error.
  const AlleleCount& target_sample_allele_count =
      allele_counts_per_sample.at(target_sample_);
  if (!nucleus::AreCanonicalBases(target_sample_allele_count.ref_base())) {
    // We don't emit calls at any site in the genome that isn't one of the
    // canonical DNA bases (one of A, C, G, or T).
    return std::nullopt;
  }

  SelectAltAllelesResult output_options = SelectAltAlleles(
      {
      .allele_counts_by_sample = allele_counts_per_sample,
      .create_complex_alleles = options_.create_complex_alleles(),
      .prev_deletion_end = prev_deletion_end
    });
  std::string ref_bases = output_options.ref_bases;
  bool complex_variant_created = output_options.complex_variant_created;

  // Determine if site is a methylated reference site based on read support
  // Exclude methylated reference sites in X or Y chromosomes from being
  // candidates for methylation-aware phasing.
  bool has_methylation = false;
  bool ref_only_site = false;
  int total_reads = target_sample_allele_count.ref_supporting_read_count();
  std::string chrom = target_sample_allele_count.position().reference_name();
  if (output_options.alt_alleles.empty() &&
      !IsExcludedMethylationContig(chrom)) {
    ref_only_site = true;
    int methyl_reads = 0;
    for (const auto& read_entry : target_sample_allele_count.read_alleles()) {
      const Allele& allele = read_entry.second;
      // Check if the allele is a reference base
      if (allele.is_methylated()) {
        methyl_reads++;
      }
    }

    // Total reads is 0 for ref sites that are unmethylated
    // Check for total_reads > 0 to avoid division by 0
    if (total_reads > 0 &&
        static_cast<double>(methyl_reads) / total_reads > 0) {
      has_methylation = true;
    }
  }
  if (output_options.alt_alleles.empty() && !KeepReferenceSite()
      && !has_methylation) {
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
  // ref_bases are calculated by SelectAltAlleles only for complex variants.
  // Otherwise it is calculated here.
  if (ref_bases.empty()) {
    ref_bases =
        CalcRefBases(target_sample_allele_count.ref_base(),
                     output_options.alt_alleles);
  }
  variant->set_reference_bases(ref_bases);
  variant->set_end(variant->start() + ref_bases.size());
  // Change the end position of the variant if there is a deletion allele.
  for (const auto& allele : output_options.alt_alleles) {
    if (allele.type() == AlleleType::DELETION) {
      prev_deletion_end = variant->start() + ref_bases.size();
      break;
    }
  }
  AddGenotypes(options_.sample_name(), {-1, -1}, variant);

  // Compute the map from read alleles to the alleles we'll use in our Variant.
  // Add the alternate alleles from our allele_map to the variant.
  const AlleleMap allele_map =
      BuildAlleleMap(target_sample_allele_count,
                     output_options.alt_alleles, ref_bases);

  // Skip next skip_next_count allele counts if comexple variant was processed.
  if (ref_bases.size() > 1 && allele_map.size() > 1 &&
        complex_variant_created) {
    skip_next_count = ref_bases.size() - 1;
  }

  for (const auto& elt : allele_map) {
    variant->add_alternate_bases(elt.second);
  }

  // If we don't have any alt_alleles, we are generating a reference site so
  // add in the kNoAltAllele.
  if (output_options.alt_alleles.empty()) {
    variant->add_alternate_bases(kNoAltAllele);
  }
  std::sort(variant->mutable_alternate_bases()->pointer_begin(),
            variant->mutable_alternate_bases()->pointer_end(),
            StringPtrLessThan());

  AddReadDepths(target_sample_allele_count, allele_map, variant);
  if (output_options.complex_variant_created) {
    AddSupportingReads(output_options.allele_counts_mod, allele_map,
                      target_sample_, &call);
  } else {
    AddSupportingReads(allele_counts_per_sample, allele_map,
                      target_sample_, &call);
  }
  if (options_.small_model_vaf_context_window_size() > 0) {
    AddAdjacentAlleleFractionsAtPosition(
        options_.small_model_vaf_context_window_size(),
        *target_sample_allele_count_iterator,
        &call);
  }

  // Add methylation information ( methylation fraction (MF) and methylation
  // depth (MD) ) to the INFO field of each call.
  ComputeMethylationStats(target_sample_allele_count, allele_map, variant);

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
    const AlleleMap& allele_map, absl::string_view target_sample,
    DeepVariantCall* call) const {
  // Iterate over each read in the allele_count, and add its name to the
  // supporting reads of for the Variant allele it supports.
  const std::string unknown_allele = kSupportingUncalledAllele;
  absl::flat_hash_map<std::string, absl::flat_hash_set<std::string>>
      alt_allele_support;
  absl::flat_hash_set<std::string> ref_support;
  for (const auto& allele_counts_entry : allele_counts) {
    const std::string& sample_name = allele_counts_entry.first;
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
        read_info->set_is_reverse_strand(allele.is_reverse_strand());
        read_info->set_sample_name(sample_name);
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
        read_info->set_is_reverse_strand(allele.is_reverse_strand());
        read_info->set_sample_name(sample_name);
      }
    }
  }
}

void VariantCaller::AddAdjacentAlleleFractionsAtPosition(
    const int window_size,
    std::vector<AlleleCount>::const_iterator
        target_sample_allele_count_iterator,
    DeepVariantCall* call) const {
  const std::vector<AlleleCount>& target_sample_allele_counts =
      allele_counters_per_sample_.at(target_sample_)->Counts();
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

// Helper function to combine the methylated reference sites and keep only the
// positive strands.
// For 5mC methylation, Pacbio marks only forward positions,
// and the reverse is assumed.
void VariantCaller::MergeMethylatedAlleleCounts(
    const std::unordered_map<std::string, AlleleCounter*>& allele_counters)
    const {
  for (const auto& sample_entry : allele_counters) {
    AlleleCounter* allele_counter = sample_entry.second;
    auto& allele_counts = allele_counter->MutableCounts();

    // Start at 1 since we are looking for the C allele preceding the G allele
    // at position i-1
    for (size_t i = 1; i < allele_counts.size(); ++i) {
      auto& allele_count = allele_counts[i];
      char ref_base = allele_count.ref_base()[0];

      // Process only G reference sites.
      if (ref_base == 'G' && IsReferenceSite(allele_count)) {
        auto& prev_allele_count = allele_counts[i - 1];

        // Ensure the preceding allele is a C and is at the genomic position
        // preceding the current G allele.
        if (prev_allele_count.ref_base()[0] != 'C' &&
            prev_allele_count.position().position() !=
                allele_count.position().position() - 1) continue;

        // Track reads that were methylated in G.
        std::vector<std::string> methylated_read_keys;

        // Remove methylation from G but store affected reads.
        for (auto& [read_key, allele] : *allele_count.mutable_read_alleles()) {
          if (allele.is_methylated()) {
            methylated_read_keys.push_back(read_key);  // Store the read key
            allele.set_is_methylated(false);  // Remove methylation
          }
        }

        // Transfer methylation to the same read keys in the C site.
        for (const std::string& read_key : methylated_read_keys) {
          auto it = prev_allele_count.mutable_read_alleles()->find(read_key);
          if (it != prev_allele_count.mutable_read_alleles()->end()) {
            it->second.set_is_methylated(true);  // Set methylation on same read
          }
        }
      }
    }
  }
}

bool VariantCaller::IsReferenceSite(const AlleleCount& allele_count) const {
  const std::string& ref_base = allele_count.ref_base();

  for (const auto& sample_entry : allele_count.sample_alleles()) {
    for (const auto& allele : sample_entry.second.alleles()) {
      if (allele.bases() != ref_base) {
        return false;
      }
    }
  }
  return true;
}

// Helper function to check if the chromosome is X or Y to exclude from
// methylation-aware phasing. Remove chr prefix if present.
bool VariantCaller::IsExcludedMethylationContig(const std::string& chrom)
    const {
  // Check if this chromosome is in the excluded list
  return std::find(options_.exclude_contigs_for_methylation_phasing().begin(),
                   options_.exclude_contigs_for_methylation_phasing().end(),
                   chrom) !=
         options_.exclude_contigs_for_methylation_phasing().end();
}

// Helper function to compute both methylation fraction and methylation depth
// Computes the methylation statistics for reference and alternate alleles.
//
// This function calculates:
// - Methylation Fraction (MF): The fraction of methylated reads for each
//   allele.
// - Methylation Depth (MD): The count of methylated reads for each allele.
//
// Low-quality reads (as determined by `allele.is_low_quality()`) are
// excluded from both calculations.
//
// @param allele_count The AlleleCount object containing read-level information.
// @param allele_map   A mapping from observed alleles to their corresponding
//                     alternate bases.
//
// @return A pair of vectors:
// - The first MF vector contains the methylation fractions for each allele.
// - The second MD vector contains the methylation depths for each allele.
//
// The vectors are guaranteed to be the same size, with each element
// corresponding to a unique allele.
// Computes and adds methylation fraction (MF) and methylation depth (MD)
// to the INFO field of the variant for a single target sample.
void VariantCaller::ComputeMethylationStats(
    const AlleleCount& target_sample_allele_count, const AlleleMap& allele_map,
    Variant* variant) const {
  // Ensure variant is not null
  CHECK(variant != nullptr);

  std::vector<double> mf_values;
  std::vector<int> md_values;

  // Compute statistics for REF allele
  int ref_methylated_count = 0, ref_total_count = 0;
  for (const auto& [read_name, allele] :
       target_sample_allele_count.read_alleles()) {
    if (allele.type() == AlleleType::REFERENCE && !allele.is_low_quality()) {
      ++ref_total_count;
      if (allele.is_methylated()) {
        ++ref_methylated_count;
      }
    }
  }

  mf_values.push_back(
      ref_total_count > 0
          ? static_cast<double>(ref_methylated_count) / ref_total_count
          : 0.0);
  md_values.push_back(ref_methylated_count);

  // Compute statistics for ALT alleles
  for (const auto& [alt_allele, alt_bases] : allele_map) {
    int alt_methylated_count = 0, alt_total_count = 0;

    for (const auto& [read_name, read_allele] :
         target_sample_allele_count.read_alleles()) {
      auto it = FindAllele(read_allele, allele_map);
      if (it != allele_map.end() && it->second == alt_bases &&
          !read_allele.is_low_quality()) {
        ++alt_total_count;
        if (read_allele.is_methylated()) {
          ++alt_methylated_count;
        }
      }
    }

    mf_values.push_back(
        alt_total_count > 0
            ? static_cast<double>(alt_methylated_count) / alt_total_count
            : 0.0);
    md_values.push_back(alt_methylated_count);
  }

  // Only add MF and MD field if at least one value is > 0
  bool has_nonzero_mf = std::any_of(
      mf_values.begin(), mf_values.end(), [](double mf) { return mf > 0.0; });

  if (has_nonzero_mf) {
    nucleus::SetInfoField(kMFFormatField, mf_values, variant->mutable_calls(0));
    nucleus::SetInfoField(kMDFormatField, md_values, variant->mutable_calls(0));
  }
}

}  // namespace multi_sample
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
