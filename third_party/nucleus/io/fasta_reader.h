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
#ifndef THIRD_PARTY_NUCLEUS_IO_FASTA_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_FASTA_READER_H_

#include <vector>
#include <unordered_map>

#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/vendor/statusor.h"

namespace nucleus {


// An FASTA reader backed by in-memory ReferenceSequence protos.
//
// FASTA files store information about DNA/RNA/Amino Acid sequences:
//
// https://en.wikipedia.org/wiki/FASTA_format
//
//
// An InMemoryRefReader provides the same API as GenomeReferenceFAI but doesn't
// fetch its data from an on-disk FASTA file but rather fetches the bases from
// an in-memory cache containing ReferenceSequence protos.
//
// In particular the GetBases(Range(chrom, start, end)) operation fetches bases
// from the tuple where chrom == chromosome, and then from the bases where the
// first base of bases starts at start. If start > 0, then the bases string is
// assumed to contain bases starting from that position in the region. For
// example, the record ('1', 10, 'ACGT') implies that
// GetBases(ranges.make_range('1', 11, 12)) will return the base 'C', as the 'A'
// base is at position 10. This makes it straightforward to cache a small region
// of a full chromosome without having to store the entire chromosome sequence
// in memory (potentially big!).
class InMemoryGenomeReference : public GenomeReference {
 public:
  // Creates a new InMemoryGenomeReference backed by ReferenceSequence protos.
  static StatusOr<std::unique_ptr<InMemoryGenomeReference>> Create(
      const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
      const std::vector<nucleus::genomics::v1::ReferenceSequence>& seqs);

  // Disable copy and assignment operations
  InMemoryGenomeReference(const InMemoryGenomeReference& other) = delete;
  InMemoryGenomeReference& operator=(const InMemoryGenomeReference&) = delete;

  const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const override {
    return contigs_;
  }

  const std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence>&
      ReferenceSequences() const {
    return seqs_;
  }

  StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const override;

 private:
  // Must use one of the static factory methods.
  explicit InMemoryGenomeReference(
      const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
      std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence>
          seqs)
      : contigs_(contigs), seqs_(seqs) {}

  const std::vector<nucleus::genomics::v1::ContigInfo> contigs_;
  const std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence>
      seqs_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_FASTA_READER_H_
