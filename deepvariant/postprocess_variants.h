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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_POSTPROCESS_VARIANTS_H_
#define LEARNING_GENOMICS_DEEPVARIANT_POSTPROCESS_VARIANTS_H_

#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

constexpr char kPhaseSetContigTag[] = "PS_CONTIG";
constexpr char kAltPhaseSetTag[] = "ALT_PS";
constexpr char kFirstVariantInPhaseSetTag[] = "FIRST_VARIANT_IN_BLOCK";
constexpr char kPhaseSetTag[] = "PS";
constexpr int kNullPhaseSetId = -1;

enum class PhaseSetStitchingStatus {
  MATCH = 0,
  SWITCH = 1,
  NOT_ENOUGH_OVERLAP = 2
};

struct VariantPhaseInformation {
  std::string phase_set_shard_id;
  std::string phase_set_region_id;
  PhaseSetStitchingStatus phase_set_stitching_status;
  bool is_first_variant_in_phase_set;
  int64_t first_variant_in_phase_set_start_position;
  bool was_phased_successfully = false;

  bool is_null() const {
    return phase_set_shard_id == std::to_string(kNullPhaseSetId) &&
           phase_set_region_id == std::to_string(kNullPhaseSetId);
  }
};

// Reads TFRecord of CallVariantsOutput protos, sort them based
// on the mapping of chromosome names to positions in FASTA in `contigs`,
// and then outputs the sorted TFRecord of CallVariantsOutput protos to
// `output_tfrecord_path`.
std::uint64_t ProcessSingleSiteCallTfRecords(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<std::string>& tfrecord_paths,
    const string& output_tfrecord_path,
    const std::vector<nucleus::genomics::v1::Range>& ranges);

// If the variant is phased, and the phase set stitching status is SWITCH,
// and the genotypes are not the same, swap the genotypes.
void MaybeSwapPhase(nucleus::genomics::v1::Variant* variant,
                    const VariantPhaseInformation& variant_phase_info);

// Reads TFRecord of Variant protos, stitch phase sets by maybe swapping the
// genotypes, and setting the PS field, which indicates the contiguous phase
// set. The output Variants are written to the output TFRecord files. The number
// of input files must match the number of output files.
void StitchPhaseSets(const std::vector<std::string>& tfrecord_paths,
                     const std::string& switches_output_path,
                     const std::vector<std::string>& output_tfrecord_paths);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_POSTPROCESS_VARIANTS_H_
