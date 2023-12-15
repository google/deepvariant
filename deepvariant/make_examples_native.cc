/*
 * Copyright 2023 Google LLC.
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

#include "deepvariant/make_examples_native.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/pileup_image_native.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Variant = nucleus::genomics::v1::Variant;
using Range = nucleus::genomics::v1::Range;

ExamplesGenerator::ExamplesGenerator(const MakeExamplesOptions& options,
                                     bool test_mode)
    : options_(options) {
  half_width_ = (options_.pic_options().width() - 1) / 2;
  if (test_mode) {
    return;
  }
  // Initialize reference reader.
  // We should always have reference name set up, except for some unit tests.
  // This code will fail if reference is not set up or cannot be loaded.
  nucleus::genomics::v1::FastaReaderOptions fasta_reader_options;
  std::string fasta_path = options_.reference_filename();
  fasta_reader_options.set_keep_true_case(false);
  std::string fai_path = absl::StrCat(fasta_path, ".fai");
  ref_reader_ = std::move(
      nucleus::IndexedFastaReader::FromFile(fasta_path, fai_path).ValueOrDie());
}


std::vector<std::vector<std::string>> ExamplesGenerator::AltAlleleCombinations(
    const Variant& variant) const {
  std::vector<std::vector<std::string>> alt_combinations;
  switch (options_.pic_options().multi_allelic_mode()) {
    case deepvariant::PileupImageOptions::UNSPECIFIED:
      LOG(FATAL) << "multi_allelic_mode cannot be UNSPECIFIED";
      break;
    case deepvariant::PileupImageOptions::NO_HET_ALT_IMAGES:
      for (const std::string& alt : variant.alternate_bases()) {
        alt_combinations.push_back({alt});
      }
      break;
    case deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES:
      {
        std::vector<std::string> alts = {variant.reference_bases()};
        alts.insert(alts.end(), variant.alternate_bases().begin(),
                    variant.alternate_bases().end());
        for (int i = 0; i < alts.size(); ++i) {
          for (int j = i + 1; j < alts.size(); ++j) {
            std::vector<std::string> one_combination;
            // Ref allele is not used in combinations.
            if (i > 0) {
              one_combination.push_back(alts[i]);
            }
            one_combination.push_back(alts[j]);
            alt_combinations.push_back(std::move(one_combination));
          }
        }
      }
      break;
    default:
      LOG(FATAL) << "Unknown value is specified for PileupImageOptions";
  }
  return alt_combinations;
}


std::string ExamplesGenerator::CreateHaplotype(const Variant& variant,
                                               const std::string& alt,
                                               int64_t* ref_start_out,
                                               int64_t* ref_end_out) const {
  int64_t var_start = variant.start();
  absl::string_view ref_bases = variant.reference_bases();
  std::string contig = variant.reference_name();
  int64_t var_end = var_start + ref_bases.size();

  std::string prefix = "";
  int64_t ref_start = std::max(var_start - half_width_, 0L);
  if (ref_start < var_start) {
    prefix =
        ref_reader_->GetBases(
          nucleus::MakeRange(contig, ref_start, var_start)).ValueOrDie();
  }
  std::string suffix = "";
  int64_t ref_end =
      std::min(ref_reader_->Contig(contig).ValueOrDie()->n_bases(),
               var_end + half_width_);
  if (ref_end > var_end) {
    suffix = ref_reader_->GetBases(
        nucleus::MakeRange(contig, var_end, ref_end)).ValueOrDie();
  }
  *ref_start_out = ref_start;
  *ref_end_out = ref_end;

  return absl::StrCat(prefix, alt, suffix);
}


// TODO There is an ongoing work to re-implement pileup channels.
// As a result this code will be obsolete once new channels code is submitted.
// This functionality is ported from deepvariant/python/clif_converters.cc.
int FillPileupArray(const std::vector<std::unique_ptr<ImageRow>>& image,
                     std::vector<unsigned char>* data) {
  int num_channels = 6;
  for (const std::unique_ptr<ImageRow>& row : image) {
    for (int i = 0; i < row->Width(); i++) {
      data->push_back(row->base[i]);
      data->push_back(row->base_quality[i]);
      data->push_back(row->mapping_quality[i]);
      data->push_back(row->on_positive_strand[i]);
      data->push_back(row->supports_alt[i]);
      data->push_back(row->matches_ref[i]);
      if (row->use_allele_frequency) {
        data->push_back(row->allele_frequency[i]);
        num_channels += 1;
      }
      if (row->add_hp_channel) {
        data->push_back(row->hp_value[i]);
        num_channels += 1;
      }
      if (!row->channels.empty()) {
        // Iterate over channels here and fill data...
        for (int j = 0; j < row->channels.size(); j++) {
          data->push_back(row->channel_data[j][i]);
          num_channels += 1;
        }
      }
    }  // for row->Width
  }  // for row
  return num_channels;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
