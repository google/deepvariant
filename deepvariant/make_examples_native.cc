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
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/pileup_image_native.h"
#include "absl/container/flat_hash_map.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/tfrecord_writer.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/example/example.pb.h"
#include "tensorflow/core/example/feature.pb.h"
#include "re2/re2.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Variant = nucleus::genomics::v1::Variant;
using Range = nucleus::genomics::v1::Range;
using Read = nucleus::genomics::v1::Read;

ExamplesGenerator::ExamplesGenerator(const MakeExamplesOptions& options,
                                     bool test_mode)
    : options_(options) {
  // Initialize samples.
  for (const auto& sample_options_for_one_sample : options_.sample_options()) {
    samples_[sample_options_for_one_sample.role()] =
        Sample(sample_options_for_one_sample);
  }
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

  // Initialize TFRecord writers for each sample.
  for (auto& [role, sample] : samples_) {
    if (sample.sample_options.skip_output_generation()) {
      continue;
    }
    // TFRecrd examples are always compressed as it is also in call_variants.
    sample.writer = nucleus::TFRecordWriter::New(
        GetExamplesFilename(options_, sample), "GZIP");
  }
}

ExamplesGenerator::~ExamplesGenerator() {
  for (auto& [role, sample] : samples_) {
    if (sample.sample_options.skip_output_generation()) {
      continue;
    }
    sample.writer->Close();
  }
}

std::string GetExamplesFilename(const MakeExamplesOptions& options,
                                const Sample& sample) {
  // If make_examples is run on a single sample or in training mode then suffix
  // is not added to the output file name.
  if (options.mode() == deepvariant::MakeExamplesOptions::TRAINING ||
      options.sample_options().size() == 1) {
    return options.examples_filename();
  } else {
    // Add suffix to the filename.
    std::string suffix = sample.sample_options.role();
    std::string file_name_with_suffix = options.examples_filename();
    RE2::Replace(
        &file_name_with_suffix, ".tfrecord",
        absl::StrCat("_", suffix, ".tfrecord"));
    return file_name_with_suffix;
  }
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
// Returns a number of channels.
int FillPileupArray(const std::vector<std::unique_ptr<ImageRow>>& image,
                     std::vector<unsigned char>* data) {
  int num_channels = 6;

  for (const std::unique_ptr<ImageRow>& row : image) {
    for (int i = 0; i < row->Width(); i++) {
      if (!row->channel_data.empty()) {
        // Iterate over channels here and fill data...
        for (int j = 0; j < row->channel_data.size(); j++) {
          data->push_back(row->channel_data[j][i]);
          num_channels += 1;
        }
      }
    }  // for row->Width
  }  // for row
  return num_channels;
}

// Calculates the variant type to be encoded in a TensorFlow example.
EncodedVariantType EncodedVariantType(const Variant& variant) {
  if (variant.reference_bases().size() == 1 &&
      variant.alternate_bases_size() >= 1) {
    bool all_alts_are_size_1 = true;
    for (const std::string& alt : variant.alternate_bases()) {
      all_alts_are_size_1 = all_alts_are_size_1 && alt.size() == 1;
    }
    if (all_alts_are_size_1) {
      return EncodedVariantType::kSnp;
    }
  }
  if (variant.reference_bases().size() > 1) {
    return EncodedVariantType::kIndel;
  }
  for (const std::string& alt : variant.alternate_bases()) {
    if (alt.size() > 1) {
      return EncodedVariantType::kIndel;
    }
  }
  return EncodedVariantType::kUnknown;
}

std::string ExamplesGenerator::EncodeExample(
    const std::vector<std::unique_ptr<ImageRow>>& image, const Variant& variant,
    const std::vector<std::string>& alt_combination) const {
  // TODO Once we know the number of channels in advance we should
  // allocate data with the correct size to avoid vector resizing when data is
  // filled.
  std::vector<unsigned char> data;
  std::array<int, 3> image_shape;
  image_shape[0] = image.size();                   // Number of rows.
  image_shape[1] = image[0]->Width();              // Width of the pileup.
  image_shape[2] = FillPileupArray(image, &data);  // Number of channels.

  // Encode alt allele indices.
  absl::flat_hash_map<std::string, int> alt_indices;
  int i = 0;
  for (const std::string& alt : variant.alternate_bases()) {
    alt_indices[alt] = i;
    i++;
  }
  std::vector<int> indices;
  indices.reserve(alt_combination.size());
  for (const std::string& alt : alt_combination) {
    indices.push_back(alt_indices[alt]);
  }
  std::string alt_indices_encoded;
  CallVariantsOutput::AltAlleleIndices alt_indices_proto;
  for (const auto& idx : indices) {
    alt_indices_proto.mutable_indices()->Add(idx);
  }
  alt_indices_proto.SerializeToString(&alt_indices_encoded);

  // Encode variant range.
  Range variant_range;
  std::string variant_range_encoded;
  variant_range.set_reference_name(variant.reference_name());
  variant_range.set_start(variant.start());
  variant_range.set_end(variant.end());
  variant_range.SerializeToString(&variant_range_encoded);

  // Ecode features of the example.
  tensorflow::Example example;
  std::string ecoded_variant;
  variant.SerializeToString(&ecoded_variant);
  (*example.mutable_features()->mutable_feature())["locus"]
      .mutable_bytes_list()
      ->add_value(variant_range_encoded);
  (*example.mutable_features()->mutable_feature())["variant/encoded"]
      .mutable_bytes_list()
      ->add_value(ecoded_variant);
  (*example.mutable_features()->mutable_feature())["variant_type"]
      .mutable_int64_list()
      ->add_value(static_cast<int64_t>(EncodedVariantType(variant)));
  (*example.mutable_features()->mutable_feature())["alt_allele_indices/encoded"]
      .mutable_bytes_list()
      ->add_value(alt_indices_encoded);
  (*example.mutable_features()->mutable_feature())["image/encoded"]
      .mutable_bytes_list()
      ->add_value(data.data(), data.size());
  for (auto dim : image_shape) {
    (*example.mutable_features()->mutable_feature())["image/shape"]
        .mutable_int64_list()
        ->add_value(dim);
  }
  (*example.mutable_features()->mutable_feature())["sequencing_type"]
      .mutable_int64_list()
      ->add_value(options_.pic_options().sequencing_type());
  std::string encoded_example;

  // Example is serialized to a string before it is written to a TFRecord.
  example.SerializeToString(&encoded_example);
  return encoded_example;
}

void ExamplesGenerator::CreateAndWritePileupExamples(
    const DeepVariantCall& candidate, const Sample& sample,
    const std::vector<int>& sample_order,
    const std::vector<InMemoryReader>& readers) {
  LOG(FATAL) << "CreateAndWritePileupExamples not implemented";
}

// Examples are generated using reads from all samples, the order of reads
// depends on the role.
// reads_per_sample contain reads for each sample. In a multisample mode reads
// are stacked according to the reads_per_sample order.
void ExamplesGenerator::WriteExamplesInRegion(
    const std::vector<nucleus::ConstProtoPtr<DeepVariantCall>>& candidates,
    const std::vector<std::vector<nucleus::ConstProtoPtr<Read>>>&
        reads_per_sample,
    const std::vector<int>& sample_order, const std::string& role) {
  // Load reads.
  std::vector<InMemoryReader> readers;
  // Cache reads passed from Python. The order of samples is preserved as passed
  // from the caller.
  readers.reserve(reads_per_sample.size());
  for (const auto& reads : reads_per_sample) {
    readers.push_back(InMemoryReader(InMemoryReader(reads)));
  }
  for (const auto& candidate : candidates) {
    CreateAndWritePileupExamples(*candidate.p_, samples_[role], sample_order,
                                 readers);
  }
}

InMemoryReader::InMemoryReader(
    const std::vector<nucleus::ConstProtoPtr<Read>>& reads)
    : reads_cache_(reads) {}

// Iterate reads_cache_ and return reads that overlap the range. In order to
// avoid copying Read protos raw pointers are returned. Since all the reads are
// allocated in Python we don't need to worry about the memory management.
std::vector<const Read*> InMemoryReader::Query(const Range& range) const {
  std::vector<const Read*> out;
  for (auto& read : reads_cache_) {
    if (nucleus::ReadOverlapsRegion(*read.p_, range)) {
      out.emplace_back(read.p_);
    }
  }
  return out;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
