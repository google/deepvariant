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

#ifndef THIRD_PARTY_NUCLEUS_IO_SAM_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_SAM_WRITER_H_

#include <memory>
#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

// A SAM/BAM/CRAM writer.
//
// SAM/BAM/CRAM files store information about a biological sequence and its
// corresponding quality scores.
//
// https://samtools.github.io/hts-specs/SAMv1.pdf
// https://samtools.github.io/hts-specs/CRAMv3.pdf
// This class converts nucleus.genomics.v1.SamHeader and
// nucleus.genomics.v1.Read to a file based on the file path passed in.
//
// This uses the htslib C API for writing NGS reads (BAM, SAM, SAM). For details
// of the API, see:
// https://github.com/samtools/htslib/tree/develop/htslib
//
class SamWriter {
 public:
  // Creates a new SamWriter writing to the file at |sam_path|, which is
  // opened and created if needed. Returns either a unique_ptr to the
  // SamWriter or a Status indicating why an error occurred.
  static StatusOr<std::unique_ptr<SamWriter>> ToFile(
      const string& sam_path,
      const nucleus::genomics::v1::SamHeader& sam_header);

  // Creates a new SamWriter writing to the file at |sam_path|, which is
  // opened and created if needed. |ref_path|, which points to an external
  // reference FASTA file, cannot be empty for CRAM files. If |embed_ref|, the
  // CRAM output file will embed the references in the output file. Returns
  // either a unique_ptr to the SamWriter or a Status indicating why an error
  // occurred.
  static StatusOr<std::unique_ptr<SamWriter>> ToFile(
      const string& sam_path, const string& ref_path, bool embed_ref,
      const nucleus::genomics::v1::SamHeader& sam_header);

  ~SamWriter();

  // Disable copy and assignment operations.
  SamWriter(const SamWriter& other) = delete;
  SamWriter& operator=(const SamWriter&) = delete;

  // Write a Read to the  file.
  // Returns Status::OK() if the write was successful; otherwise the status
  // provides information about what error occurred.
  ::nucleus::Status Write(const nucleus::genomics::v1::Read& read);
  ::nucleus::Status WritePython(
      const ConstProtoPtr<const nucleus::genomics::v1::Read>& wrapped) {
    return Write(*(wrapped.p_));
  }

  // Close the underlying resource descriptors. Returns Status::OK() if the
  // close was successful; otherwise the status provides information about what
  // error occurred.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support. Do
  // not use it!
  void PythonEnter() const {}

 private:
  class NativeHeader;
  class NativeFile;
  class NativeBody;
  // Private constructor; use ToFile to safely create a SamWriter.
  SamWriter(std::unique_ptr<NativeFile> file,
            std::unique_ptr<NativeHeader> header);

  // A pointer to the htslib file used to access the SAM/BAM/CRAM data.
  std::unique_ptr<NativeFile> native_file_;

  // A htslib header data structure obtained by parsing the header of this file.
  std::unique_ptr<NativeHeader> native_header_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_SAM_WRITER_H_
