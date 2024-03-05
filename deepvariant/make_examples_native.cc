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

#include "deepvariant/alt_aligned_pileup_lib.h"
#include "deepvariant/pileup_image_native.h"
#include "absl/container/flat_hash_map.h"
#include "absl/log/check.h"
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
    : options_(options), pileup_image_(options_.pic_options()) {
  // Initialize samples.
  for (const auto& sample_options_for_one_sample : options_.sample_options()) {
    samples_[sample_options_for_one_sample.role()] =
        Sample(sample_options_for_one_sample);
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
                                               absl::string_view alt,
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

// Converting a vector<ImageRow> into a 3-d vector where channels are stored in
// the lowest dimension. Special processing is required for alt aligned
// channels.
// * --alt_aligned_pileup=diff_channels Add channel 5 (base differs
// from ref) of both alts as channels
// * --alt_aligned_pileup=base_channels Add channel 0 (bases ATCG) of both alts
// as channels.
// * --alt_aligned_pileup=rows alt aligned data is stored as extra rows.
//
void FillPileupArray(
    const std::vector<std::unique_ptr<ImageRow>>& image,
    const std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_image,
    const AltAlignedPileup alt_aligned_representation,
    std::vector<unsigned char>* pileup_array) {
  // TODO Although channels 0 and 5 should not change we need to set
  // channel number programmatically.
  int alt_channel_index = 0;
  switch (alt_aligned_representation) {
    case AltAlignedPileup::kDiffChannels:
      alt_channel_index = 5;
      break;
    case AltAlignedPileup::kBaseChannels:
      alt_channel_index = 0;
      break;
    case AltAlignedPileup::kRows:
      LOG(FATAL) << "Alt alignment as row is not implemented.";
    default:
      alt_channel_index = 0;
  }
  for (int row = 0; row < image.size(); row++) {
    for (int column = 0; column < image[row]->Width(); column++) {
      if (!image[row]->channel_data.empty()) {
        // Lower dimension is a channel data. Here we iterate all channels to
        // fill one position of the pileup image.
        for (int channel = 0; channel < image[row]->channel_data.size();
             channel++) {
          pileup_array->push_back(image[row]->channel_data[channel][column]);
        }
        // Fill alt aligned channels if needed.
        if (alt_aligned_representation == AltAlignedPileup::kBaseChannels ||
            alt_aligned_representation == AltAlignedPileup::kDiffChannels) {
          CHECK_EQ(alt_image.size(), 2);
          unsigned char alt_aligned_channel_1_value = 0;
          if (!alt_image[0].empty() &&
              !alt_image[0][row]->channel_data.empty()) {
            alt_aligned_channel_1_value =
                alt_image[0][row]->channel_data[alt_channel_index][column];
          }
          // Fill alt aligned channel 1.
          pileup_array->push_back(alt_aligned_channel_1_value);
          // Fill alt aligned channel 2.
          if (alt_image[1].empty() || alt_image[1][row]->channel_data.empty()) {
            // Fill with alt aligned channel 1 if alt2 is empty.
            pileup_array->push_back(alt_aligned_channel_1_value);
          } else {
            pileup_array->push_back(
                alt_image[1][row]->channel_data[alt_channel_index][column]);
          }
        }  // if need_alt_alignment
      }
    }  // for row->Width
  }  // for row
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
    const std::vector<std::unique_ptr<ImageRow>>& image,
    const std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_image,
    const Variant& variant,
    const std::vector<std::string>& alt_combination) const {
  // TODO Once we know the number of channels in advance we should
  // allocate data with the correct size to avoid vector resizing when data is
  // filled.
  std::vector<unsigned char> data;
  std::array<int, 3> image_shape;
  FillPileupArray(image, alt_image, alt_aligned_pileup_, &data);
  image_shape[0] = image.size();                   // Number of rows.
  image_shape[1] = image[0]->Width();              // Width of the pileup.
  image_shape[2] = options_.pic_options().channels().size();  // Num channels.

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

namespace {
  bool HasAtLeastOneNonSingleBaseAllele(const Variant& variant) {
    for (const std::string& alt : variant.alternate_bases()) {
      if (alt.size() > 1) return true;
    }
    return false;
  }
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
// ths haplotype instead of the reference. BuildPileupForOneSample is called
// the same way as for a normal pileup which means no special alt aligned
// processing is required in pileup_image_native library.
void ExamplesGenerator::CreateAltAlignedImages(
    const DeepVariantCall& candidate,
    const std::vector<std::string>& alt_combination,
    const std::vector<Read>& trimmed_reads, int sample_order,
    const Range& region,
    std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_images) {
  const auto& variant = candidate.variant();
  int image_start_pos = variant.start() - half_width_;
  int alt_image_num = 0;
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
    for (const auto& read : realigned_reads) {
      realigned_reads_ptrs.push_back(&read);
    }
    // Build alt pileup with all channels. Although, we only need one
    // channel it is easier to call BuildPileupForOneSample for all
    // channels. The runtime penalty is very small (~1.5% of runtime).
    auto alt_image = pileup_image_.BuildPileupForOneSample(
        candidate, haplotype.substr(0, options_.pic_options().width()),
        realigned_reads_ptrs, image_start_pos, alt_combination,
        options_.sample_options(sample_order));
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
    const std::vector<int>& sample_order,
    const std::vector<InMemoryReader>& readers) {
  const auto& variant = candidate.variant();
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
  bool need_alt_alignment = NeedAltAlignment(variant);
  for (const std::vector<std::string>& alt_combination :
        AltAlleleCombinations(variant)) {
    std::vector<std::unique_ptr<ImageRow>> ref_images;
    std::vector<std::vector<std::unique_ptr<ImageRow>>> alt_images(2);
    for (int this_sample_order : sample_order) {
      // All reads are trimmed if trim_reads_for_pileup is set. This allows for
      // a better runtime performance when long reads data is used.
      std::vector<Read> trimmed_reads;
      std::vector<const Read*> trimmed_reads_ptrs;
      if (options_.trim_reads_for_pileup()) {
        trimmed_reads = TrimReads(
            readers[this_sample_order].Query(region),
            CalculateAlignmentRegion(variant, half_width_, *ref_reader_));
        trimmed_reads_ptrs.reserve(trimmed_reads.size());
        for (const auto& read : trimmed_reads) {
          trimmed_reads_ptrs.push_back(&read);
        }
      }
      // If trim_reads_for_pileup is set then trimmed_reads are used.
      auto ref_image = pileup_image_.BuildPileupForOneSample(
          candidate, reference_bases,
          options_.trim_reads_for_pileup()
              ? trimmed_reads_ptrs
              : readers[this_sample_order].Query(region),
          image_start_pos, alt_combination,
          options_.sample_options(this_sample_order));
      // Collect rows from all samples in ref_images.
      for (auto& row : ref_image) {
        ref_images.push_back(std::move(row));
      }
      if (need_alt_alignment) {
        CreateAltAlignedImages(candidate, alt_combination, trimmed_reads,
                               this_sample_order, region, alt_images);
      }
    }
    sample.writer->WriteRecord(
        EncodeExample(ref_images, alt_images, variant, alt_combination));
  }
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
  // Cache reads passed from Python. The order of samples is preserved as
  // passed from the caller.
  readers.reserve(reads_per_sample.size());
  for (const auto& reads : reads_per_sample) {
    readers.push_back(InMemoryReader(InMemoryReader(reads)));
  }
  for (const auto& candidate : candidates) {
    CreateAndWriteExamplesForCandidate(*candidate.p_, samples_[role],
                                        sample_order, readers);
  }
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

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
