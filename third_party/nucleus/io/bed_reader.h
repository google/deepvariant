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

#ifndef THIRD_PARTY_NUCLEUS_IO_BED_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_BED_READER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/bed.pb.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "third_party/nucleus/core/status.h"

namespace nucleus {


// Alias for the abstract base class for BED record iterables.
using BedIterable = Iterable<nucleus::genomics::v1::BedRecord>;

// A BED reader.
//
// BED files are flexible stores of information about a genome annotation track.
//
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
//
// This class provides a method to iterate through a BED file.
//
// The objects returned by iterate() are nucleus.genomics.v1.BedRecord
// objects parsed from the BED records in the file.
//
// Note: Only tab-delimited BED files are supported, for ease of future support
// for tabix-indexed BED file querying.
//
class BedReader : public Reader {
 public:
  // Creates a new BedReader reading reads from the BED file bed_path.
  //
  // bed_path must point to an existing BED formatted file.
  //
  // Returns a StatusOr that is OK if the BedReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<BedReader>> FromFile(
      const string& bed_path,
      const nucleus::genomics::v1::BedReaderOptions& options);

  ~BedReader();

  // Disable copy and assignment operations.
  BedReader(const BedReader& other) = delete;
  BedReader& operator=(const BedReader&) = delete;

  // Gets all of the BED records in this file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<BedIterable>> Iterate() const;

  // Close the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.
  void PythonEnter() const {}

  // Get the options controlling the behavior of this BedReader.
  const nucleus::genomics::v1::BedReaderOptions& Options() const {
    return options_;
  }

  // Returns the header that tracks the number of fields in each record in the
  // reader.
  const nucleus::genomics::v1::BedHeader& Header() const { return header_; }

  // Returns OK if the input numTokens equals num_fields in the header.
  ::nucleus::Status Validate(const int numTokens) const;

 private:
  // Private constructor; use FromFile to safely create a BedReader from a
  // file.
  BedReader(std::unique_ptr<TextReader> text_reader,
            const nucleus::genomics::v1::BedReaderOptions& options,
            const nucleus::genomics::v1::BedHeader& header);

  // Our options that control the behavior of this class.
  const nucleus::genomics::v1::BedReaderOptions options_;

  // The header that tracks the number of fields in each record in the file.
  const nucleus::genomics::v1::BedHeader header_;

  // A pointer to a raw TextReader object.
  std::unique_ptr<TextReader> text_reader_;

  // Allow BedIterable objects to access fp_.
  friend class BedFullFileIterable;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_BED_READER_H_
