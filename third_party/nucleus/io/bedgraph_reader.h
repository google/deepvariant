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

#ifndef THIRD_PARTY_NUCLEUS_IO_BEDGRAPH_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_BEDGRAPH_READER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/bedgraph.pb.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/vendor/statusor.h"

namespace nucleus {

// Alias for the abstract base class for record iterables.
using BedGraphIterable = Iterable<nucleus::genomics::v1::BedGraphRecord>;

// A BedGraph reader.
//
// BedGraph files stores data values associated with genome sequences in a track
// format.
//
// https://genome.ucsc.edu/goldenpath/help/bedgraph.html
//
// This class provides a method to iterate through a BedGraph file.
//
// The objects returned by iterate() are nucleus.genomics.v1.BedGraphRecord
// objects parsed from the BedGraph records in the file.
//
class BedGraphReader : public Reader {
 public:
  // Creates a new BedGraphReader reading reads from the BedGraph file at
  // |bedgraph_path|.
  //
  // Returns a StatusOr that is OK if the BedGraphReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<BedGraphReader>> FromFile(
      const string& bedgraph_path);

  ~BedGraphReader();

  // Disable copy and assignment operations.
  BedGraphReader(const BedGraphReader& other) = delete;
  BedGraphReader& operator=(const BedGraphReader&) = delete;

  // Gets all of the BedGraph records in this file in order. Returns an OK
  // status if the iterable can be constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<BedGraphIterable>> Iterate() const;

  // Closes the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.
  void PythonEnter() const {}

 private:
  // Private constructor. Use FromFile to safely create a BedGraphReader from a
  // file.
  BedGraphReader(std::unique_ptr<TextReader> text_reader);

  // A pointer to a raw TextReader object.
  std::unique_ptr<TextReader> text_reader_;

  // Allow BedGraphIterable objects to access fp_.
  friend class BedGraphFullFileIterable;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_BEDGRAPH_READER_H_
