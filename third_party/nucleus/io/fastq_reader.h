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
 */

#ifndef THIRD_PARTY_NUCLEUS_IO_FASTQ_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_FASTQ_READER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/fastq.pb.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/vendor/statusor.h"

namespace nucleus {

// Alias for the abstract base class for FASTQ record iterables.
using FastqIterable = Iterable<nucleus::genomics::v1::FastqRecord>;

// A FASTQ reader.
//
// FASTQ files store information about a biological sequence and its
// corresponding quality scores.
//
// https://en.wikipedia.org/wiki/FASTQ_format
//
// This class provides a method to iterate through a FASTQ file.
//
// The objects returned by iterate() are nucleus.genomics.v1.FastqRecord
// objects parsed from the FASTQ records in the file.
//
class FastqReader : public Reader {
 public:
  // Creates a new FastqReader reading reads from the FASTQ file fastq_path.
  //
  // fastq_path must point to an existing FASTQ formatted file.
  //
  // Returns a StatusOr that is OK if the FastqReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<FastqReader>> FromFile(
      const string& fastq_path,
      const nucleus::genomics::v1::FastqReaderOptions& options);

  ~FastqReader();

  // Disable copy and assignment operations.
  FastqReader(const FastqReader& other) = delete;
  FastqReader& operator=(const FastqReader&) = delete;

  // Gets all of the FASTQ records in this file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<FastqIterable>> Iterate() const;

  // Close the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it! Returns a Status indicating whether the enter was successful.
  ::nucleus::Status PythonEnter() const { return ::nucleus::Status(); }

  // Get the options controlling the behavior of this FastqReader.
  const nucleus::genomics::v1::FastqReaderOptions& Options() const {
    return options_;
  }

 private:
  // Private constructor; use FromFile to safely create a FastqReader from a
  // file.
  FastqReader(std::unique_ptr<TextReader> text_reader,
              const nucleus::genomics::v1::FastqReaderOptions& options);

  // Populates the four string  pointers with values from the input file.
  ::nucleus::Status Next(string* header, string* sequence, string* pad,
                         string* quality) const;

  // Our options that control the behavior of this class.
  const nucleus::genomics::v1::FastqReaderOptions options_;

  // Underlying file reader.
  std::unique_ptr<TextReader> text_reader_;

  // Give Iterator classes access to Next().
  friend class FastqFullFileIterable;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_FASTQ_READER_H_
