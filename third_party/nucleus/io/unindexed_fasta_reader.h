/*
 * Copyright 2018 Google LLC.
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

#ifndef THIRD_PARTY_NUCLEUS_IO_UNINDEXED_FASTA_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_UNINDEXED_FASTA_READER_H_

#include <memory>
#include <string>
#include <vector>

#include "absl/types/optional.h"
#include "htslib/faidx.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

// A FASTA reader that is not backed by a htslib FAI index.
//
// FASTA files store information about DNA/RNA/Amino Acid sequences:
//
// https://en.wikipedia.org/wiki/FASTA_format
//
// This reader is for FASTA files that contain many small records and are
// explicitly not indexed. The FASTA files can be optionally block-gzipped
// compressed.
//
// This class provides methods to iterate through a the FASTA records but
// doesn't support query() for the bases spanning a specific region on the
// genome.
//
// The (name, bases) tuple returned by iterate() are strings containing the
// bases in uppercase.
class UnindexedFastaReader : public GenomeReference {
 public:
  // Creates a new GenomeReference backed by the FASTA file fasta_path.
  //
  // Returns this newly allocated UnindexedFastaReader object, passing ownership
  // to the caller via a unique_ptr.
  static StatusOr<std::unique_ptr<UnindexedFastaReader>> FromFile(
      const string& fasta_path);

  ~UnindexedFastaReader();

  // Disable copy and assignment operations
  UnindexedFastaReader(const UnindexedFastaReader& other) = delete;
  UnindexedFastaReader& operator=(const UnindexedFastaReader&) = delete;

  const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const override;

  StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const override;

  StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>> Iterate()
      const override;

  // Close the underlying resource descriptors.
  tensorflow::Status Close() override;

 private:
  // Allow iteration to access the underlying reader.
  friend class UnindexedFastaReaderIterable;

  // Must use one of the static factory methods.
  UnindexedFastaReader(std::unique_ptr<TextReader> text_reader);

  const std::vector<nucleus::genomics::v1::ContigInfo> contigs_;

  // Underlying file reader.
  std::unique_ptr<TextReader> text_reader_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_UNINDEXED_FASTA_READER_H_
