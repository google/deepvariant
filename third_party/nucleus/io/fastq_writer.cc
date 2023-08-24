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

// Implementation of fastq_writer.h
#include "third_party/nucleus/io/fastq_writer.h"

#include <utility>

#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/fastq.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"

namespace nucleus {

// 256 KB write buffer.
constexpr int WRITER_BUFFER_SIZE = 256 * 1024;

// -----------------------------------------------------------------------------
//
// Writer for FASTQ formats containing NGS reads.
//
// -----------------------------------------------------------------------------

StatusOr<std::unique_ptr<FastqWriter>> FastqWriter::ToFile(
    const string& fastq_path,
    const nucleus::genomics::v1::FastqWriterOptions& options) {
  StatusOr<std::unique_ptr<TextWriter>> text_writer =
      TextWriter::ToFile(fastq_path);
  NUCLEUS_RETURN_IF_ERROR(text_writer.status());
  return absl::WrapUnique(
      new FastqWriter(text_writer.ConsumeValueOrDie(), options));
}

FastqWriter::FastqWriter(
    std::unique_ptr<TextWriter> text_writer,
    const nucleus::genomics::v1::FastqWriterOptions& options)
    : options_(options), text_writer_(std::move(text_writer)) {}

FastqWriter::~FastqWriter() {
  if (text_writer_) {
    NUCLEUS_CHECK_OK(Close());
  }
}

::nucleus::Status FastqWriter::Close() {
  if (!text_writer_)
    return ::nucleus::FailedPrecondition(
        "Cannot close an already closed FastqWriter");
  // Close the file pointer we have been writing to.
  ::nucleus::Status close_status = text_writer_->Close();
  text_writer_ = nullptr;
  return close_status;
}

::nucleus::Status FastqWriter::Write(
    const nucleus::genomics::v1::FastqRecord& record) {
  if (!text_writer_)
    return ::nucleus::FailedPrecondition(
        "Cannot write to closed FASTQ stream.");
  string out = "@";
  absl::StrAppend(&out, record.id());
  if (!record.description().empty()) {
    absl::StrAppend(&out, " ", record.description());
  }
  absl::StrAppend(&out, "\n", record.sequence(), "\n+\n", record.quality(),
                  "\n");
  NUCLEUS_RETURN_IF_ERROR(text_writer_->Write(out));

  return ::nucleus::Status();
}

}  // namespace nucleus
