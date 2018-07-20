/*
 * Copyright 2018 Google Inc.
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
 *
 */

#include "third_party/nucleus/io/fasta_reader.h"

#include <algorithm>

#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/errors.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::ReferenceSequence;

// Initializes an InMemoryGenomeReference from contigs and seqs.
//
// contigs is a vector describing the "contigs" of this GenomeReference. These
// should include only the contigs present in seqs. A ContigInfo object for a
// contig `chrom` should describe the entire chromosome `chrom` even if the
// corresponding ReferenceSequence only contains a subset of the bases.
//
// seqs is a vector where each element describes a region of the genome we are
// caching in memory and will use to provide bases in the query() operation.
//
// Note that only a single ReferenceSequence for each contig is currently
// supported.
//
// There should be exactly one ContigInfo for each reference_name referred to
// across all ReferenceSequences, and no extra ContigInfos.
StatusOr<std::unique_ptr<InMemoryGenomeReference>>
InMemoryGenomeReference::Create(
      const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
      const std::vector<nucleus::genomics::v1::ReferenceSequence>& seqs) {
  std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence> seqs_map;

  for (const auto& seq : seqs) {
    if (seq.region().reference_name().empty() || seq.region().start() < 0 ||
        seq.region().start() > seq.region().end()) {
      return tensorflow::errors::InvalidArgument(
          "Malformed region ", seq.region().ShortDebugString());
    }

    const size_t region_len = seq.region().end() - seq.region().start();
    if (region_len != seq.bases().length()) {
      return tensorflow::errors::InvalidArgument(
          "Region size = ", region_len, " not equal to bases.length() ",
          seq.bases().length());
    }

    auto insert_result = seqs_map.emplace(seq.region().reference_name(), seq);
    if (!insert_result.second) {
      return tensorflow::errors::InvalidArgument(
          "Each ReferenceSequence must be on a different chromosome but "
          "multiple ones were found on ", seq.region().reference_name());
    }
  }

  return std::unique_ptr<InMemoryGenomeReference>(
      new InMemoryGenomeReference(contigs, seqs_map));
}

StatusOr<string> InMemoryGenomeReference::GetBases(const Range& range) const {
  if (!IsValidInterval(range))
    return tensorflow::errors::InvalidArgument("Invalid interval: ",
                                               range.ShortDebugString());

  const ReferenceSequence& seq = seqs_.at(range.reference_name());

  if (range.start() < seq.region().start() ||
      range.end() > seq.region().end()) {
    return tensorflow::errors::InvalidArgument(
        "Cannot query range=", range.ShortDebugString(),
        " as this InMemoryRefReader only has bases in the interval=",
        seq.region().ShortDebugString());
  }
  const int64 pos = range.start() - seq.region().start();
  const int64 len = range.end() - range.start();
  return seq.bases().substr(pos, len);
}

}  // namespace nucleus
