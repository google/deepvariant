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

#ifndef THIRD_PARTY_NUCLEUS_IO_GFF_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_GFF_READER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

const nucleus::genomics::v1::GffReaderOptions kDefaultGffReaderOptions{};

// Alias for the abstract base class for GFF record iterables.
using GffIterable = Iterable<nucleus::genomics::v1::GffRecord>;

class GffReader : public Reader {
 public:
  // Creates a new GffReader reading reads from the GFF file gff_path.
  //
  // gff_path must point to an existing GFF formatted file (or gzipped
  // equivalent).
  //
  // The GFF format is described here:
  // https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
  //
  // Returns a StatusOr that is OK if the GffReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<GffReader>> FromFile(
      const string& gff_path,
      const nucleus::genomics::v1::GffReaderOptions& options =
          kDefaultGffReaderOptions);

  ~GffReader() = default;

  // Disable copy and assignment operations.
  GffReader(const GffReader& other) = delete;
  GffReader& operator=(const GffReader&) = delete;

  // Gets all of the GFF records in this file in order.
  // Returns an OK status if the iterable can be constructed, or not
  // OK otherwise.  Iteration is over proto records of type
  // nucleus.genomics.v1.GffRecord
  StatusOr<std::shared_ptr<GffIterable>> Iterate() const;

  // Closes the underlying resource descriptors. Returns a Status to
  // indicate if everything went OK with the close.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.
  void PythonEnter() const {}

  // Get the options controlling the behavior of this GffReader.
  const nucleus::genomics::v1::GffReaderOptions& Options() const {
    return options_;
  }

  // Returns the header that tracks the number of fields in each record in the
  // reader.
  const nucleus::genomics::v1::GffHeader& Header() const { return header_; }

 private:
  // Private constructor used by FromFile factory.
  GffReader(std::unique_ptr<TextReader> text_reader,
            const nucleus::genomics::v1::GffReaderOptions& options,
            const nucleus::genomics::v1::GffHeader& header);

  // A pointer to a raw TextReader object.
  std::unique_ptr<TextReader> text_reader_;

  // Options controlling the behavior of this class.
  const nucleus::genomics::v1::GffReaderOptions options_;

  // The GFF header, reflecting how to interpret fields.
  const nucleus::genomics::v1::GffHeader header_;

  // Allow iteration to access the underlying reader.
  friend class GffFullFileIterable;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_GFF_READER_H_
