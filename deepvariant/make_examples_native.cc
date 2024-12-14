/*
 * Copyright 2024 Google LLC.
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
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "deepvariant/alt_aligned_pileup_lib.h"
#include "deepvariant/pileup_image_native.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/stream_examples.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
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

constexpr int kDefaultMinimumReadOverlap = 15;

// ExamplesGenerator is inialized once per each make_examples instance.
ExamplesGenerator::ExamplesGenerator(
    const MakeExamplesOptions& options,
    const std::unordered_map<std::string, std::string>& example_filenames,
    bool test_mode)
    : options_(options), pileup_image_(options_.pic_options()) {
  // Initialize samples.
  for (const auto& sample_options_for_one_sample : options_.sample_options()) {
    samples_[sample_options_for_one_sample.role()] =
        Sample(sample_options_for_one_sample);
    samples_[sample_options_for_one_sample.role()].needs_alt_alignment =
        SampleNeedsAltAlignment(sample_options_for_one_sample);
    // Blank out specified channels.
    for (int i = 0;
         i < sample_options_for_one_sample.channels_enum_to_blank_size(); ++i) {
      samples_[sample_options_for_one_sample.role()]
          .channels_enum_to_blank.insert(
              sample_options_for_one_sample.channels_enum_to_blank(i));
    }
    // Blank out all channels for the specified variant types.
    for (int i = 0;
         i < sample_options_for_one_sample.variant_types_to_blank_size(); ++i) {
      switch (sample_options_for_one_sample.variant_types_to_blank(i)) {
        case SampleOptions::VARIANT_TYPE_SNP:
          samples_[sample_options_for_one_sample.role()]
              .encoded_variant_types_to_blank.insert(EncodedVariantType::kSnp);
          break;
        case SampleOptions::VARIANT_TYPE_INDEL:
          samples_[sample_options_for_one_sample.role()]
              .encoded_variant_types_to_blank.insert(
                  EncodedVariantType::kIndel);
          break;
        default:
          LOG(WARNING) << "Unknown variant type to blank will be ignored: "
                       << sample_options_for_one_sample.variant_types_to_blank(
                              i);
          break;
      }
    }
  }
  half_width_ = (options_.pic_options().width() - 1) / 2;
  absl::flat_hash_map<absl::string_view, AltAlignedPileup>
      alt_aligned_pileup_map({
          {"none", AltAlignedPileup::kNone},
          {"base_channels", AltAlignedPileup::kBaseChannels},
          {"diff_channels", AltAlignedPileup::kDiffChannels},
          {"rows", AltAlignedPileup::kRows},
      });
  auto it =
      alt_aligned_pileup_map.find(options_.pic_options().alt_aligned_pileup());
  if (it == alt_aligned_pileup_map.end()) {
    LOG(FATAL) << "Unknown value is specified for alt_aligned_pileup";
  } else {
    alt_aligned_pileup_ = it->second;
  }
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


  if (options_.stream_examples()) {
    stream_examples_ = std::make_unique<StreamExamples>(options_,
                                                        alt_aligned_pileup_);
  } else {
    // Initialize TFRecord writers for each sample.
    // TFRecord writers are not used if examples are streamed.
    for (auto& [role, sample] : samples_) {
      if (sample.sample_options.skip_output_generation()) {
        continue;
      }
      auto it = example_filenames.find(role);
      if (it == example_filenames.end()) {
        LOG(INFO) << "Example filename not found for role: " << role;
        continue;
      }
      // TFRecrd examples are always compressed as it is also in
      // call_variants.
      sample.writer = nucleus::TFRecordWriter::New(it->second, "GZIP");
      if (sample.writer == nullptr) {
        LOG(FATAL) << "Failed to create TFRecord writer for " << it->second;
      }
    }
  }
}

ExamplesGenerator::~ExamplesGenerator() {
  if (!options_.stream_examples()) {
    for (auto& [role, sample] : samples_) {
      if (sample.sample_options.skip_output_generation()) {
        continue;
      }
      if (sample.writer != nullptr) {
        sample.writer->Close();
      }
    }
  }
}

std::string GetExamplesFilename(const MakeExamplesOptions& options,
                                const Sample& sample, bool add_suffix) {
  // If make_examples is run on a single sample or in training mode then suffix
  // is not added to the output file name.
  if (options.mode() == deepvariant::MakeExamplesOptions::TRAINING ||
      !add_suffix) {
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
                                               absl::string_view alt,
                                               int64_t* ref_start_out,
                                               int64_t* ref_end_out) const {
  int64_t var_start = variant.start();
  absl::string_view ref_bases = variant.reference_bases();
  const std::string& contig = variant.reference_name();
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

namespace {
bool HasAtLeastOneNonSingleBaseAllele(const Variant& variant) {
  for (const std::string& alt : variant.alternate_bases()) {
    if (alt.size() > 1) return true;
  }
  return false;
}

void UpdateStats(enum EncodedVariantType variant_type,
                 const VariantLabel* label, int label_value,
                 std::unordered_map<std::string, int>& stats) {
  stats["n_examples"] += 1;
  if (variant_type == EncodedVariantType::kIndel) {
    stats["n_indels"] += 1;
  } else {
    stats["n_snps"] += 1;
  }

  if (label != nullptr) {
    stats["n_class_0"] += (label_value == 0 ? 1 : 0);
    stats["n_class_1"] += (label_value == 1 ? 1 : 0);
    stats["n_class_2"] += (label_value == 2 ? 1 : 0);
    stats["n_non_denovo"] += (label->is_denovo ? 0 : 1);
    stats["n_denovo"] += (label->is_denovo ? 1 : 0);
  }
}

std::string EncodeAltAlleles(const Variant& variant,
                             absl::Span<const std::string> alt_combination,
                             absl::flat_hash_set<int>* alt_indices_set) {
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
    if (alt_indices_set != nullptr) {
      alt_indices_set->insert(alt_indices[alt]);
    }
  }
  CallVariantsOutput::AltAlleleIndices alt_indices_proto;
  for (const auto& idx : indices) {
    alt_indices_proto.mutable_indices()->Add(idx);
  }
  std::string alt_indices_encoded;
  alt_indices_proto.SerializeToString(&alt_indices_encoded);
  return alt_indices_encoded;
}

std::string EncodeVariant(const Variant& variant, const VariantLabel* label) {
  std::string encoded_variant;
  if (label == nullptr) {
    variant.SerializeToString(&encoded_variant);
  } else {  // For train examples the variant from label is used.
    label->variant.SerializeToString(&encoded_variant);
  }
  return encoded_variant;
}

}  // namespace

std::string ExamplesGenerator::EncodeExample(
    absl::Span<const std::unique_ptr<ImageRow>> image,
    absl::Span<const std::vector<std::unique_ptr<ImageRow>>> alt_image,
    const Variant& variant, absl::Span<const std::string> alt_combination,
    std::unordered_map<std::string, int>& stats, std::vector<int>& image_shape,
    const std::unique_ptr<VariantLabel>& label) const {
  // if AltAlignedPileup::kRows is set then number of
  // rows equals: (image.size + alt_image_1.size + alt_image_2.size) or
  //   image.size * 3.
  if (alt_aligned_pileup_ == AltAlignedPileup::kRows) {
    // Width of pileup for AltAlignedPileup::kRows
    image_shape[0] = image.size() * 3;
  } else {
    image_shape[0] = image.size();  // Number of rows.
  }
  image_shape[1] = image[0]->Width();              // Width of the pileup.
  image_shape[2] = options_.pic_options().channels().size();  // Num channels.

  std::vector<unsigned char> data(
      image_shape[0] * image_shape[1] * image_shape[2], 0);
  FillPileupArray(image, alt_image, alt_aligned_pileup_, &data, data.size());

  // Encode alt allele indices.
  absl::flat_hash_set<int> alt_indices_set;
  std::string alt_indices_encoded =
      EncodeAltAlleles(variant, alt_combination, &alt_indices_set);

  // Encode variant range.
  Range variant_range;
  std::string variant_range_encoded;
  std::ostringstream s;
  // The string literal form looks like:
  //   reference_name:start+1-end
  // since start and end are zero-based inclusive (start) and exclusive (end),
  // while the literal form is one-based inclusive on both ends.
  s << variant.reference_name() << ":" << variant.start() + 1 << "-"
    << variant.end();
  variant_range_encoded.assign(s.str());

  // Encode features of the example.
  tensorflow::Example example;
  enum EncodedVariantType variant_type = EncodedVariantType(variant);
  (*example.mutable_features()->mutable_feature())["locus"]
      .mutable_bytes_list()
      ->add_value(std::move(variant_range_encoded));
  (*example.mutable_features()->mutable_feature())["variant/encoded"]
      .mutable_bytes_list()
      ->add_value(EncodeVariant(variant, label.get()));
  (*example.mutable_features()->mutable_feature())["variant_type"]
      .mutable_int64_list()
      ->add_value(static_cast<int64_t>(variant_type));
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

  // Set the label if it is provided.
  int label_value = 0;
  if (label != nullptr) {
    label_value = label->LabelForAltAlleles(alt_indices_set);
    (*example.mutable_features()->mutable_feature())["label"]
        .mutable_int64_list()
        ->add_value(label_value);
  }

  // Set de novo feature if de novo regions are provided.
  if (!options_.denovo_regions_filename().empty() && label != nullptr) {
    (*example.mutable_features()->mutable_feature())["denovo_label"]
        .mutable_int64_list()
        ->add_value(label->is_denovo);
  }

  // Example is serialized to a string before it is written to a TFRecord.
  std::string encoded_example;
  example.SerializeToString(&encoded_example);
  UpdateStats(variant_type, label.get(), label_value, stats);
  return encoded_example;
}

bool ExamplesGenerator::SampleNeedsAltAlignment(
    const SampleOptions& sample_options) const {
  // Right now, if any of the alt aligned channels are listed in
  // `channels_enum_to_blank` of the sample option, we'll skip alt alignment.
  // In the future, we may consider a more fine-grained approach.
  const std::vector<DeepVariantChannelEnum> alt_aligned_channels = {
    DeepVariantChannelEnum::CH_DIFF_CHANNELS_ALTERNATE_ALLELE_1,
    DeepVariantChannelEnum::CH_DIFF_CHANNELS_ALTERNATE_ALLELE_2,
    DeepVariantChannelEnum::CH_BASE_CHANNELS_ALTERNATE_ALLELE_1,
    DeepVariantChannelEnum::CH_BASE_CHANNELS_ALTERNATE_ALLELE_2
  };
  if (!sample_options.channels_enum_to_blank().empty()) {
    for (const auto& channel_enum : alt_aligned_channels) {
      for (const auto& channel_enum_to_blank :
               sample_options.channels_enum_to_blank()) {
        if (channel_enum_to_blank == channel_enum) {
          return false;
        }
      }
    }
  }
  return true;
}

bool ExamplesGenerator::NeedAltAlignment(const Variant& variant) const {
  if (options_.pic_options().alt_aligned_pileup() == "none") {
    return false;
  }
  if (options_.pic_options().types_to_alt_align() == "all") {
    return true;
  }
  if (options_.pic_options().types_to_alt_align() == "indels") {
    return variant.reference_bases().size() > 1 ||
            HasAtLeastOneNonSingleBaseAllele(variant);
  }
  return false;
  }

std::string ExamplesGenerator::GetReferenceBasesForPileup(
    const Variant& variant) const {
  int start = variant.start() - half_width_;
  int end = start + options_.pic_options().width();
  Range region;
  region.set_reference_name(variant.reference_name());
  region.set_start(start);
  region.set_end(end);
  return GetReferenceBases(region);
}

std::string ExamplesGenerator::GetReferenceBases(const Range& region) const {
  if (ref_reader_->IsValidInterval(region)) {
    return ref_reader_->GetBases(region).ValueOrDie();
  }
  return "";
}

// Two pileup images are built, one for each alt allele. In these pileups the
// reference is replaced with the artificial haplotype, that is created by
// concatenating ref prefix, alt allele, and ref suffix. Reads are aligned to
// this haplotype instead of the reference. BuildPileupForOneSample is called
// the same way as for a normal pileup which means no special alt aligned
// processing is required in pileup_image_native library.
void ExamplesGenerator::CreateAltAlignedImages(
    const DeepVariantCall& candidate,
    const std::vector<std::string>& alt_combination,
    const std::vector<Read>& trimmed_reads, int sample_order,
    const Range& region,
    std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_images,
    const float mean_coverage, std::vector<int64_t>* original_start_positions) {
  const auto& variant = candidate.variant();
  int image_start_pos = variant.start() - half_width_;
  int alt_image_num = 0;
  CHECK(original_start_positions == nullptr ||
        original_start_positions->size() == trimmed_reads.size());

  // It is assumed that there can be 0 - 2 alt alleles in alt_combination.
  CHECK_LE(alt_combination.size(), 2);
  for (const std::string& alt : alt_combination) {
    int64_t ref_start = 0;
    int64_t ref_end = 0;
    // Create haplotype by combining reference and alt bases in the
    // middle.
    std::string haplotype = CreateHaplotype(variant, alt, &ref_start, &ref_end);
    if (haplotype.size() < options_.pic_options().width()) {
      // We need haplotype length to be the size of pilup width. In some
      // cases, for example, when variant is too close to the beginning
      // of a contig the haplotype length would be smaller than pilup
      // width. In this case we do not fill alt aligned channels.
      break;
    }
    // Realign reads to the haplotype.
    std::vector<Read> realigned_reads = RealignReadsToHaplotype(
        haplotype, trimmed_reads, region.reference_name(), ref_start, ref_end,
        *(this->ref_reader_), this->options_);
    std::vector<const Read*> realigned_reads_ptrs;
    realigned_reads_ptrs.reserve(realigned_reads.size());

    // Some reads cannot be aligned to haplotype. In that case empty read is
    // inserted in realigned_reads vector. We need to remove those and modify
    // original_start_positions accordingly.
    int index = 0;
    std::vector<int64_t> start_positions_after_realignment;
    for (const auto& read : realigned_reads) {
      if (!read.aligned_sequence().empty()) {
        realigned_reads_ptrs.push_back(&read);
        start_positions_after_realignment.push_back(
            original_start_positions->at(index));
      }
      index++;
    }
    // Build alt pileup with all channels. Although, we only need one
    // channel it is easier to call BuildPileupForOneSample for all
    // channels. The runtime penalty is very small (~1.5% of runtime).
    auto sample_it =
        samples_.find(options_.sample_options(sample_order).role());
    // Get the channels to blank for this sample if any.
    auto blank_channels_set = sample_it->second.channels_enum_to_blank;
    if (sample_it->second.encoded_variant_types_to_blank.contains(
            EncodedVariantType(variant))) {
      auto channel_enum_descriptor = DeepVariantChannelEnum_descriptor();
      for (int i = 0; i < channel_enum_descriptor->value_count(); ++i) {
        blank_channels_set.insert(static_cast<DeepVariantChannelEnum>(
            channel_enum_descriptor->value(i)->number()));
      }
    }
    auto alt_image = pileup_image_.BuildPileupForOneSample(
        candidate, haplotype.substr(0, options_.pic_options().width()),
        realigned_reads_ptrs, image_start_pos, alt_combination,
        options_.sample_options(sample_order), mean_coverage,
        &start_positions_after_realignment, blank_channels_set);
    // move alt_image to alt_images[2] array.
    for (auto& row : alt_image) {
      alt_images[alt_image_num].push_back(std::move(row));
    }
    alt_image_num++;
  }  // for (alt in alt_combination)
}

// This function is tested by integration test in make_examples_test using the
// real data.
// For multisample cases image is created for each sample separately and all
// images are stacked into one.
void ExamplesGenerator::CreateAndWriteExamplesForCandidate(
    const DeepVariantCall& candidate, const Sample& sample,
    absl::Span<const int> sample_order,
    absl::Span<const InMemoryReader> readers,
    std::unordered_map<std::string, int>& stats, std::vector<int>& image_shape,
    absl::Span<const float> mean_coverage_per_sample,
    const std::unique_ptr<VariantLabel>& label) {
  const auto& variant = candidate.variant();
  const auto encoded_variant_type = EncodedVariantType(variant);
  int image_start_pos = variant.start() - half_width_;
  // Pileup range.
  Range region;
  region.set_reference_name(candidate.variant().reference_name());
  region.set_start(candidate.variant().start() -
                    options_.pic_options().read_overlap_buffer_bp());
  region.set_end(candidate.variant().end() +
                  options_.pic_options().read_overlap_buffer_bp());
  std::string reference_bases = GetReferenceBasesForPileup(variant);
  if (reference_bases.empty()) {
    // We are at the edge of the contig, example cannot be created.
    return;
  }
  bool variant_needs_alt_alignment = NeedAltAlignment(variant);
  bool use_trimmed_reads =
      options_.trim_reads_for_pileup() || variant_needs_alt_alignment ||
      sample.sample_options.keep_only_window_spanning_reads();
  int min_trimming_overlap = kDefaultMinimumReadOverlap;
  for (const std::vector<std::string>& alt_combination :
        AltAlleleCombinations(variant)) {
    std::vector<std::unique_ptr<ImageRow>> ref_images;
    std::vector<std::vector<std::unique_ptr<ImageRow>>> alt_images(2);
    for (int this_sample_order : sample_order) {
      // Implementing the logic from internal#comment5.
      CHECK(this_sample_order < options_.sample_options().size());
      if (options_.sample_options(this_sample_order)
              .keep_only_window_spanning_reads()) {
        min_trimming_overlap = options_.pic_options().width();
      } else {
        min_trimming_overlap = kDefaultMinimumReadOverlap;
      }
      // All reads are trimmed if trim_reads_for_pileup is set. This allows for
      // a better runtime performance when long reads data is used.
      std::vector<Read> trimmed_reads;
      std::vector<const Read*> trimmed_reads_ptrs;
      std::vector<int64_t> original_start_positions = {};
      // We always trim reads if alt alignment is needed.
      if (use_trimmed_reads) {
        trimmed_reads = TrimReads(
            /*reads=*/readers[this_sample_order].Query(region),
            /*region=*/
            CalculateAlignmentRegion(variant, half_width_, *ref_reader_),
            /*original_alignment_positions=*/original_start_positions,
            /*min_overlap=*/min_trimming_overlap);
        trimmed_reads_ptrs.reserve(trimmed_reads.size());
        for (const auto& read : trimmed_reads) {
          trimmed_reads_ptrs.push_back(&read);
        }
      }
      auto mean_coverage = mean_coverage_per_sample[this_sample_order];
      // If trim_reads_for_pileup is set then trimmed_reads are used.
      auto sample_it =
          samples_.find(options_.sample_options(this_sample_order).role());
      // Get the channels to blank for this sample if any.
      auto blank_channels_set = sample_it->second.channels_enum_to_blank;
      if (sample_it->second.encoded_variant_types_to_blank.contains(
              EncodedVariantType(variant))) {
        auto channel_enum_descriptor = DeepVariantChannelEnum_descriptor();
        for (int i = 0; i < channel_enum_descriptor->value_count(); ++i) {
          blank_channels_set.insert(static_cast<DeepVariantChannelEnum>(
              channel_enum_descriptor->value(i)->number()));
        }
      }
      auto ref_image = pileup_image_.BuildPileupForOneSample(
          candidate, reference_bases,
          use_trimmed_reads ? std::move(trimmed_reads_ptrs)
                            : readers[this_sample_order].Query(region),
          image_start_pos, alt_combination,
          options_.sample_options(this_sample_order), mean_coverage,
          &original_start_positions, blank_channels_set);
      // Collect rows from all samples in ref_images.
      for (auto& row : ref_image) {
        ref_images.push_back(std::move(row));
      }
      bool sample_needs_alt_alignment = sample_it->second.needs_alt_alignment;
      if (variant_needs_alt_alignment) {
        if (sample_needs_alt_alignment) {
          CreateAltAlignedImages(candidate, alt_combination, trimmed_reads,
                                 this_sample_order, region, alt_images,
                                 mean_coverage, &original_start_positions);
        } else {
          // If the sample does not need alt alignment, we still want to fill
          // the corresponding area with 0s.
          // NOTE: The easiest way to simulate this behavior is to use the same
          // function, but to give it an empty list of reads.
          // This might not be the cleanest way to write the logic, but it is
          // a lot easier for human to understand.
          CreateAltAlignedImages(candidate, alt_combination, {},
                                 this_sample_order, region, alt_images,
                                 mean_coverage, nullptr);
        }
      }
    }
    if (options_.stream_examples()) {
      stream_examples_->StreamExample(
          ref_images, alt_images, alt_aligned_pileup_,
          EncodeAltAlleles(variant, alt_combination, nullptr),
          EncodeVariant(variant, nullptr));
      UpdateStats(encoded_variant_type, nullptr, 0, stats);
    } else {
      sample.writer->WriteRecord(EncodeExample(ref_images, alt_images, variant,
                                               alt_combination, stats,
                                               image_shape, label));
    }
  }
}

// Examples are generated using reads from all samples, the order of reads
// depends on the role.
// reads_per_sample contain reads for each sample. In a multisample mode reads
// are stacked according to the reads_per_sample order.
std::unordered_map<std::string, int> ExamplesGenerator::WriteExamplesInRegion(
    absl::Span<const nucleus::ConstProtoPtr<DeepVariantCall>> candidates,
    absl::Span<const std::vector<nucleus::ConstProtoPtr<Read>>>
        reads_per_sample,
    const std::vector<int>& sample_order, const std::string& role,
    // const std::vector<VariantLabel>& labels,
    const std::vector<float>& mean_coverage_per_sample,
    std::vector<int>* image_shape) {
  CHECK(labels_.empty() || candidates.size() == labels_.size());
  auto sample_it = samples_.find(role);
  CHECK(sample_it != samples_.end()) << "Role " << role << " not found.";
  CHECK(sample_it->second.writer != nullptr || options_.stream_examples())
      << "Role " << role << " does not have a writer.";

  // image_shape is the return parameter that is passed to Python. The memory
  // is handled by Python.
  image_shape->resize(3);
  // Load reads.
  std::vector<InMemoryReader> readers;
  std::unordered_map<std::string, int> stats;
  // Cache reads passed from Python. The order of samples is preserved as
  // passed from the caller.
  readers.reserve(reads_per_sample.size());
  for (const auto& reads : reads_per_sample) {
    readers.push_back(InMemoryReader(InMemoryReader(reads)));
  }
  // Acquire mutex until all examples are written or no space left in the
  // buffer.
  if (options_.stream_examples()) {
    stream_examples_->StartStreaming();
  }

  for (int i = 0; i < candidates.size(); i++) {
    if (labels_.empty()) {
      CreateAndWriteExamplesForCandidate(
          *(candidates[i].p_), sample_it->second, sample_order, readers, stats,
          *image_shape, mean_coverage_per_sample, nullptr);
    } else {
      CreateAndWriteExamplesForCandidate(
          *(candidates[i].p_), sample_it->second, sample_order, readers, stats,
          *image_shape, mean_coverage_per_sample, labels_[i]);
    }
  }
  // Write zero to indicate the end of the examples only if we had any examples
  // created.
  if (options_.stream_examples()) {
    stream_examples_->EndStreaming(!candidates.empty());
  }

  labels_.clear();
  return stats;
}

InMemoryReader::InMemoryReader(
    const std::vector<nucleus::ConstProtoPtr<Read>>& reads)
    : reads_cache_(reads) {}

// Iterate reads_cache_ and return reads that overlap the range. In order to
// avoid copying Read protos raw pointers are returned. Since all the reads
// are allocated in Python we don't need to worry about the memory management.
std::vector<const Read*> InMemoryReader::Query(const Range& range) const {
  std::vector<const Read*> out;
  for (auto& read : reads_cache_) {
    if (nucleus::ReadOverlapsRegion(*read.p_, range)) {
      out.emplace_back(read.p_);
    }
  }
  return out;
}

// Implementation of VariantLabel. This functionality is moved from
// variant_labeler.py.
int VariantLabel::LabelForAltAlleles(
    const absl::flat_hash_set<int>& alt_indices_set) const {
  int label_value = 0;

  for (int i = 0; i < genotype.size(); i++) {
    if (genotype[i] == 0) {
      continue;
    }
    if (alt_indices_set.contains(genotype[i] - 1)) {
      label_value++;
    }
  }
  CHECK(label_value >= 0 && label_value <= 2);
  return label_value;
}

absl::string_view CustomizedClassesLabel::GetClassStatus(
    const ::google::protobuf::Map<std::string, ::nucleus::genomics::v1::ListValue>&
        info_field) const {
  if (info_field.find(info_field_name) == info_field.end()) {
    LOG(FATAL) << "Cannot create class labels: VCF file does not contain INFO/"
               << info_field_name << " field";
  }

  absl::string_view class_status =
      info_field.at(info_field_name).values(0).string_value();

  if (classes_dict.find(std::string(class_status)) == classes_dict.end()) {
    std::stringstream all_classes;
    for (const auto& [class_status, label] : classes_dict) {
      all_classes << class_status << ", ";
    }
    LOG(FATAL) << "class_status status unknown: " << class_status
               << " Known status: " << all_classes.str();
  }
  return class_status;
}

// Implementation of CustomizedClassesLabel. The code is adopted from
// customized_class_labeler.py.
// THe function works the following way:
// 1. If the truth variant is {0,0} then the label is 0.
// 2. If the ref bases of the candidate and the truth doesn't match, return
//    0 (ref).
// 3. class status of truth is returned if first alt allele of truth is in
//    variant's alt alleles set. For example, if first alt allele of truth is
//    'AC' and variant's alt alleles are ['T', 'AC'] then there is a match and
//    the class from truth_variant is returned.
int CustomizedClassesLabel::LabelForAltAlleles(
    const absl::flat_hash_set<int>& alt_indices_set) const {
  int label_value = 0;
  // If truth variant is {0,0} then the label is 0.
  if (!truth_variant.calls().empty() &&
      truth_variant.calls()[0].genotype()[0] == 0 &&
      truth_variant.calls()[0].genotype()[1] == 0) {
    return 0;
  }

  // If the ref of the candidate and the truth doesn't match, return 0 (ref).
  if (truth_variant.reference_bases() != variant.reference_bases()) {
    return 0;
  }

  absl::string_view true_class_status = GetClassStatus(truth_variant.info());
  auto truth_alt = truth_variant.alternate_bases()[0];
  // Default is label 0. Usually reference.
  // Note that this logic below might not be the best when
  // `alt_alleles_indices` is a composite one, like [0, 1]. For now we'll
  // return the corresponding label if any of them matches truth_alt.
  for (int ind : alt_indices_set) {
    if (variant.alternate_bases()[ind] == truth_alt) {
      // allele in called variant is the same as truth_alt
      label_value = classes_dict.at(std::string(true_class_status));
    }
  }

  return label_value;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
