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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_TESTING_UTILS_H_
#define LEARNING_GENOMICS_DEEPVARIANT_TESTING_UTILS_H_

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#endif  // LEARNING_GENOMICS_DEEPVARIANT_TESTING_UTILS_H_

#include <memory>
#include <string>
#include <vector>

#include "deepvariant/pileup_image_native.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Default sample name for unit testing.
inline constexpr absl::string_view kSampleName = "MySampleName";
// Default chromosome for unit testing.
inline constexpr absl::string_view kChr = "chr1";
// Default start of a variant used in unit tests.
inline constexpr int64_t kStart = 10;

// Creates reference contig from a given sequence. This is useful to create
// InMemoryReferenceReader for testing purposes.
void CreateTestSeq(const std::string& name, int pos_in_fasta,
                   int range_start, int range_end,
                   const std::string& bases,
                   std::vector<nucleus::genomics::v1::ContigInfo>* contigs,
                   std::vector<nucleus::genomics::v1::ReferenceSequence>* seqs);

// Creates Read proto.
nucleus::genomics::v1::Read MakeRead(
    absl::string_view chromosome, int start, const std::string& bases,
    const std::vector<std::string>& cigar_elements,
    absl::string_view read_name,
    int hp_tag = -1);

// Creates a simple variant for unit testing.
nucleus::genomics::v1::Variant MakeVariant(
    absl::string_view ref, absl::Span<const absl::string_view> alts,
    int64_t start = kStart);

DeepVariantCall MakeDeepVariantCall(
    const nucleus::genomics::v1::Variant& variant);

std::unique_ptr<ImageRow> MakeImageRow(
    const std::vector<std::vector<unsigned char>>& data, int width,
    int num_channels);

PileupImageOptions MakeDefaultPileupImageOptions(int width,
                                                 int height,
                                                 int ref_band_height);
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
