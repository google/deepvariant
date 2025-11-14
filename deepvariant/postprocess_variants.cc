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

#include "deepvariant/postprocess_variants.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/io/compression.h"
#include "tensorflow/core/lib/io/record_reader.h"
#include "tensorflow/core/lib/io/record_writer.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

namespace {

void SortSingleSiteCalls(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    std::vector<CallVariantsOutput>* calls) {
  std::vector<CallVariantsOutput> output;
  if (calls->empty()) {
    return;
  }
  //   Create the mapping from from contig to pos_in_fasta.
  std::map<std::string, int> contig_name_to_pos_in_fasta =
      nucleus::MapContigNameToPosInFasta(contigs);
  std::stable_sort(calls->begin(), calls->end(),
            [&contig_name_to_pos_in_fasta](const CallVariantsOutput& a,
                                           const CallVariantsOutput& b) {
              return nucleus::CompareVariants(a.variant(), b.variant(),
                                              contig_name_to_pos_in_fasta);
            });
}

// Loads phasing info from merge_reads output.
// Returns a map from (phase_set_shard_id, phase_set_region_id) to a
// PhaseSetStitchingStatus integer.
std::map<std::pair<std::string, std::string>, int> LoadPhasingInfo(
    const std::string& switches_output_path) {
  std::map<std::pair<std::string, std::string>, int> phasing_info;
  std::ifstream infile(switches_output_path);
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    std::vector<std::string> parts = absl::StrSplit(line, '\t');
    CHECK_EQ(parts.size(), 3) << "Invalid line in switches file: " << line;
    int switch_status;
    CHECK(absl::SimpleAtoi(parts[2], &switch_status))
        << "Invalid switch status: " << parts[2];
    phasing_info[{parts[0], parts[1]}] = switch_status;
  }
  return phasing_info;
}

// Looks up the FIRST_VARIANT_IN_BLOCK info field in the variant.
bool IsFirstVariantInPhaseSet(const Variant& variant) {
  if (variant.info().find(kFirstVariantInPhaseSetTag) == variant.info().end()) {
    return false;
  }
  return variant.info().at(kFirstVariantInPhaseSetTag).values(0).bool_value();
}

// Parses the PS_CONTIG info field in the variant.
bool GetShardAndRegionFromPsContig(
    const Variant& variant, std::pair<std::string, std::string>* phase_set) {
  if (variant.info().find(kPhaseSetContigTag) != variant.info().end() &&
      variant.info().at(kPhaseSetContigTag).values_size() > 0) {
    std::vector<std::string> parts = absl::StrSplit(
        variant.info().at(kPhaseSetContigTag).values(0).string_value(), '-');
    CHECK_EQ(parts.size(), 2)
        << "Invalid PS_CONTIG: "
        << variant.info().at(kPhaseSetContigTag).values(0).string_value();
    *phase_set = {parts[0], parts[1]};
    return true;
  }
  return false;
}

// Returns the phase set information for the given variant.
VariantPhaseInformation GetVariantPhaseInformation(
    const Variant& variant,
    const std::map<std::pair<std::string, std::string>, int>&
        stitching_status_by_phase_set,
    const VariantPhaseInformation& previous_variant_phase_info) {
  std::pair<std::string, std::string> phase_set;
  bool has_ps_contig = GetShardAndRegionFromPsContig(variant, &phase_set);
  if (!has_ps_contig) {
    return previous_variant_phase_info;
  }

  std::string phase_set_shard_id = phase_set.first;
  std::string phase_set_region_id = phase_set.second;
  int stitching_status_integer = 0;
  if (stitching_status_by_phase_set.find(phase_set) !=
      stitching_status_by_phase_set.end()) {
    stitching_status_integer = stitching_status_by_phase_set.at(phase_set);
  }
  PhaseSetStitchingStatus stitching_status =
      static_cast<PhaseSetStitchingStatus>(stitching_status_integer);

  if (previous_variant_phase_info.is_null()) {
    return VariantPhaseInformation{phase_set_shard_id, phase_set_region_id,
                                   stitching_status,   true,
                                   variant.start(),    false};
  }

  if (phase_set_shard_id == previous_variant_phase_info.phase_set_shard_id &&
      phase_set_region_id == previous_variant_phase_info.phase_set_region_id) {
    VariantPhaseInformation new_info = previous_variant_phase_info;
    if (previous_variant_phase_info.was_phased_successfully) {
      new_info.is_first_variant_in_phase_set = false;
    }
    return new_info;
  }

  // Start a new phase set.
  if (IsFirstVariantInPhaseSet(variant) ||
      stitching_status == PhaseSetStitchingStatus::NOT_ENOUGH_OVERLAP) {
    return VariantPhaseInformation{phase_set_shard_id, phase_set_region_id,
                                   stitching_status,   true,
                                   variant.start(),    false};
  }

  return VariantPhaseInformation{
      previous_variant_phase_info.phase_set_shard_id,
      previous_variant_phase_info.phase_set_region_id,
      stitching_status,
      false,
      previous_variant_phase_info.first_variant_in_phase_set_start_position,
      false};
}

}  // namespace

void MaybeSwapPhase(nucleus::genomics::v1::Variant* variant,
                    const VariantPhaseInformation& variant_phase_info) {
  if (!(variant->info().count(kPhaseSetContigTag) &&
        variant->info().count(kAltPhaseSetTag)) ||
      !variant->calls(0).is_phased()) {
    return;
  }

  VariantCall* call = variant->mutable_calls(0);
  if (variant_phase_info.phase_set_stitching_status ==
          PhaseSetStitchingStatus::SWITCH &&
      variant->calls(0).genotype(0) != variant->calls(0).genotype(1)) {
    // Swap the genotypes.
    std::vector<int> ordered_gt = {call->genotype(0), call->genotype(1)};
    std::swap(ordered_gt[0], ordered_gt[1]);
    call->clear_genotype();
    for (int gt : ordered_gt) {
      call->add_genotype(gt);
    }
  }
  call->set_is_phased(true);
  // Set PS info field.
  nucleus::genomics::v1::ListValue ps_info;
  nucleus::genomics::v1::Value* ps_value = ps_info.add_values();
  ps_value->set_int_value(
      variant_phase_info.first_variant_in_phase_set_start_position + 1);
  (*call->mutable_info())[kPhaseSetTag] = ps_info;
}

std::uint64_t ProcessSingleSiteCallTfRecords(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<std::string>& tfrecord_paths,
    const string& output_tfrecord_path,
    const std::vector<nucleus::genomics::v1::Range>& ranges) {
  std::vector<CallVariantsOutput> single_site_calls;
  tensorflow::Env* env = tensorflow::Env::Default();
  for (const string& tfrecord_path : tfrecord_paths) {
    std::unique_ptr<tensorflow::RandomAccessFile> read_file;
    TF_CHECK_OK(env->NewRandomAccessFile(tfrecord_path, &read_file));
    const char* const option = tensorflow::io::compression::kGzip;
    tensorflow::io::RecordReader reader(
        read_file.get(),
        tensorflow::io::RecordReaderOptions::CreateRecordReaderOptions(option));

    std::uint64_t offset = 0;
    tensorflow::tstring data;
    LOG(INFO) << "Read from: " << tfrecord_path;
    while (reader.ReadRecord(&offset, &data).ok()) {
      CallVariantsOutput single_site_call;
      QCHECK(single_site_call.ParseFromArray(data.data(), data.length()))
          << "Failed to parse CallVariantsOutput";
      // Here we assume each variant has only 1 call.
      QCHECK_EQ(single_site_call.variant().calls_size(), 1);
      if (ranges.empty() ||
          nucleus::RangesContainVariant(ranges, single_site_call.variant())) {
        single_site_calls.push_back(std::move(single_site_call));
      }
    }
    if (tfrecord_paths.size() > 1) {
      LOG(INFO) << "Done reading: " << tfrecord_path
                << ". #entries in single_site_calls = "
                << std::to_string(single_site_calls.size());
    }
  }
  LOG(INFO) << "Total #entries in single_site_calls = "
            << std::to_string(single_site_calls.size());
  VLOG(3) << "Start SortSingleSiteCalls";
  SortSingleSiteCalls(contigs, &single_site_calls);
  VLOG(3) << "Done SortSingleSiteCalls";

  // Write sorted calls to output_tfrecord_path.
  std::unique_ptr<tensorflow::WritableFile> output_file;
  TF_CHECK_OK(tensorflow::Env::Default()->NewWritableFile(output_tfrecord_path,
                                                          &output_file));
  tensorflow::io::RecordWriter output_writer(output_file.get());
  for (const auto& single_site_call : single_site_calls) {
    absl::Status writer_status =
        output_writer.WriteRecord(single_site_call.SerializeAsString());
    QCHECK(writer_status.ok())
        << "Failed to write serialized proto to output_writer. "
        << "Status = " << writer_status.message();
  }
  TF_CHECK_OK(output_writer.Flush()) << "Failed to flush the output writer.";
  return single_site_calls.size();
}

void StitchPhaseSets(absl::Span<const std::string> tfrecord_paths,
                     const std::string& switches_output_path,
                     absl::Span<const std::string> output_tfrecord_paths) {
  std::map<std::pair<std::string, std::string>, int>
      stitching_status_by_phase_set = LoadPhasingInfo(switches_output_path);
  LOG(INFO) << "Loaded " << stitching_status_by_phase_set.size()
            << " entries from switches file.";
  CHECK_EQ(tfrecord_paths.size(), output_tfrecord_paths.size());

  VariantPhaseInformation previous_variant_phase_info{
      std::to_string(kNullPhaseSetId),
      std::to_string(kNullPhaseSetId),
      PhaseSetStitchingStatus::MATCH,
      false,
      -1,
      false};

  tensorflow::Env* env = tensorflow::Env::Default();

  // The tfrecord_paths must be in sorted order.
  int num_variants = 0;
  for (int i = 0; i < tfrecord_paths.size(); ++i) {
    const string& tfrecord_path = tfrecord_paths[i];
    const string& output_tfrecord_path = output_tfrecord_paths[i];
    std::unique_ptr<tensorflow::RandomAccessFile> read_file;
    TF_CHECK_OK(env->NewRandomAccessFile(tfrecord_path, &read_file));
    const char* const option = tensorflow::io::compression::kNone;
    tensorflow::io::RecordReader reader(
        read_file.get(),
        tensorflow::io::RecordReaderOptions::CreateRecordReaderOptions(option));

    std::unique_ptr<tensorflow::WritableFile> output_file;
    TF_CHECK_OK(tensorflow::Env::Default()->NewWritableFile(
        output_tfrecord_path, &output_file));
    tensorflow::io::RecordWriter output_writer(output_file.get());

    std::uint64_t offset = 0;
    tensorflow::tstring data;
    while (reader.ReadRecord(&offset, &data).ok()) {
      Variant variant;
      QCHECK(variant.ParseFromArray(data.data(), data.length()));

      VariantPhaseInformation variant_phase_info = GetVariantPhaseInformation(
          variant, stitching_status_by_phase_set, previous_variant_phase_info);
      MaybeSwapPhase(&variant, variant_phase_info);
      previous_variant_phase_info = variant_phase_info;
      if (variant.calls(0).is_phased()) {
        previous_variant_phase_info.was_phased_successfully = true;
      }
      absl::Status writer_status =
          output_writer.WriteRecord(variant.SerializeAsString());
      QCHECK(writer_status.ok())
          << "Failed to write serialized proto to output_writer. "
          << "Status = " << writer_status.message();
      num_variants++;
    }
    TF_CHECK_OK(output_writer.Flush()) << "Failed to flush the output writer.";
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
