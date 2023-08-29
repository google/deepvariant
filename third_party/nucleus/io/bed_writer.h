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

#ifndef THIRD_PARTY_NUCLEUS_IO_BED_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_BED_WRITER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/text_writer.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/bed.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

// A BED writer, allowing us to write BED files.
//
// BED files flexibly store annotation information about a reference genome.
//
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
class BedWriter {
 public:
  // Creates a new BedWriter writing to the file at bed_path, which is
  // opened and created if needed. Returns either a unique_ptr to the
  // BedWriter or a Status indicating why an error occurred.
  static StatusOr<std::unique_ptr<BedWriter>> ToFile(
      const string& bed_path, const nucleus::genomics::v1::BedHeader& header,
      const nucleus::genomics::v1::BedWriterOptions& options);

  ~BedWriter();

  // Disable copy and assignment operations.
  BedWriter(const BedWriter& other) = delete;
  BedWriter& operator=(const BedWriter&) = delete;

  // Write a BedRecord to the BED file.
  // Returns Status::OK() if the write was successful; otherwise the status
  // provides information about what error occurred.
  ::nucleus::Status Write(const nucleus::genomics::v1::BedRecord& record);
  ::nucleus::Status WritePython(
      const ConstProtoPtr<const nucleus::genomics::v1::BedRecord>& wrapped) {
    return Write(*(wrapped.p_));
  }

  // Close the underlying resource descriptors. Returns Status::OK() if the
  // close was successful; otherwise the status provides information about what
  // error occurred.
  ::nucleus::Status Close();

  // Provide access to the header.
  const nucleus::genomics::v1::BedHeader& Header() const { return header_; }

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  void PythonEnter() const {}

 private:
  // Private constructor; use ToFile to safely create a BedWriter.
  BedWriter(std::unique_ptr<TextWriter> text_writer,
            const nucleus::genomics::v1::BedHeader& header,
            const nucleus::genomics::v1::BedWriterOptions& options);

  // The header of the BED file.
  const nucleus::genomics::v1::BedHeader header_;

  // Our options that control the behavior of this class.
  const nucleus::genomics::v1::BedWriterOptions options_;

  // Underlying file writer.
  std::unique_ptr<TextWriter> text_writer_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_BED_WRITER_H_
