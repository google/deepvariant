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

#ifndef THIRD_PARTY_NUCLEUS_IO_GFF_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_GFF_WRITER_H_

#include <memory>
#include <string>

#include "third_party/nucleus/io/text_writer.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

const nucleus::genomics::v1::GffWriterOptions kDefaultGffWriterOptions{};

// A GFF writer class, allowing us to write GFF (text) files.
// The GFF format is described here:
//   https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
class GffWriter {
 public:
  // Creates a new GffWriter writing to the file at gff_path, which is
  // opened and created if needed. Returns either a unique_ptr to the
  // GffWriter or a Status indicating why an error occurred.  If
  // gff_path ends in ".gz", the resulting file will be GZIP compressed.
  static StatusOr<std::unique_ptr<GffWriter>> ToFile(
      const string& gff_path, const nucleus::genomics::v1::GffHeader& header,
      const nucleus::genomics::v1::GffWriterOptions& options =
          kDefaultGffWriterOptions);

  ~GffWriter() = default;

  // Disable copy and assignment operations.
  GffWriter(const GffWriter& other) = delete;
  GffWriter& operator=(const GffWriter&) = delete;

  // Writes a GffRecord to the GFF file.
  // Returns Status::OK() if the write was successful; otherwise the status
  // provides information about what error occurred.
  ::nucleus::Status Write(const nucleus::genomics::v1::GffRecord& record);
  ::nucleus::Status WritePython(
      const ConstProtoPtr<const nucleus::genomics::v1::GffRecord>& wrapped) {
    return Write(*(wrapped.p_));
  }

  // Closes the underlying resource descriptors. Returns Status::OK() if the
  // close was successful; otherwise the status provides information about what
  // error occurred.
  ::nucleus::Status Close();

  // Provides access to the header.
  const nucleus::genomics::v1::GffHeader& Header() const { return header_; }

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  void PythonEnter() const {}

 private:
  // Private constructor; use ToFile to safely create a GffWriter.
  GffWriter(std::unique_ptr<TextWriter> text_writer,
            const nucleus::genomics::v1::GffHeader& header,
            const nucleus::genomics::v1::GffWriterOptions& options);

  // The header of the GFF file.
  const nucleus::genomics::v1::GffHeader header_;

  // Our options that control the behavior of this class.
  const nucleus::genomics::v1::GffWriterOptions options_;

  // Underlying file writer.
  std::unique_ptr<TextWriter> text_writer_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_GFF_WRITER_H_
