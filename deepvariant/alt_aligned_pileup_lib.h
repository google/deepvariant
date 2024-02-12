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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_ALT_ALIGNED_PILEUP_LIB_H_
#define LEARNING_GENOMICS_DEEPVARIANT_ALT_ALIGNED_PILEUP_LIB_H_

#include <cstdint>
#include <memory>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Trim long read starting from ref_start position and spanning for ref_length.
// Note, that starting position relative to the read as well as new read length
// have to be calculated.
// read_start - read relative position matching the ref_start.
// new_read_length - read_length after trimming.
void TrimCigar(
    const google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit>& cigar,
    int64_t ref_start, int64_t ref_length,
    google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit>* new_cigar,
    int64_t& read_start, int64_t& new_read_length);

// Trim long read to region.
nucleus::genomics::v1::Read TrimRead(
    const nucleus::genomics::v1::Read& read,
    const nucleus::genomics::v1::Range& region);

// Calculates the alignment region to match the width of a pileup with the
// given variant in the middle.
nucleus::genomics::v1::Range CalculateAlignmentRegion(
    const nucleus::genomics::v1::Variant& variant, int half_width,
    const nucleus::GenomeReference& ref_reader);

// Realign reads to a given haplotype.
// Returns a vector of new reads.
std::vector<nucleus::genomics::v1::Read> RealignReadsToHaplotype(
    absl::string_view haplotype,
    const std::vector<nucleus::genomics::v1::Read>& reads,
    absl::string_view contig,  // Chromosome name for the haplotype.
    int64_t ref_start,   // Start position of the haplotype relative to the ref.
    int64_t ref_end,     // End position of the haplotype relative to the ref.
    const nucleus::GenomeReference& ref_reader,
    const MakeExamplesOptions& options);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_ALT_ALIGNED_PILEUP_LIB_H_
