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
 */

// Implementation of fastq_writer.h
#include "third_party/nucleus/io/fastq_writer.h"

#include "absl/memory/memory.h"
#include "third_party/nucleus/protos/fastq.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/zlib_compression_options.h"
#include "third_party/nucleus/vendor/zlib_outputbuffer.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/file_system.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

namespace tf = tensorflow;

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
  std::unique_ptr<tensorflow::WritableFile> fp;
  if (!tf::Env::Default()->NewWritableFile(fastq_path.c_str(), &fp).ok()) {
    return tf::errors::Unknown(
        tf::strings::StrCat("Could not open fastq_path ", fastq_path));
  }
  bool isCompressed = EndsWith(fastq_path, ".gz");
  auto writer =
      absl::WrapUnique(new FastqWriter(std::move(fp), options, isCompressed));
  return std::move(writer);
}

FastqWriter::FastqWriter(
    std::unique_ptr<tensorflow::WritableFile> fp,
    const nucleus::genomics::v1::FastqWriterOptions& options,
    const bool isCompressed)
    : options_(options), raw_file_(std::move(fp)), isCompressed_(isCompressed) {
  if (isCompressed_) {
    auto zwriter = new tf::io::ZlibOutputBuffer(
        raw_file_.get(), WRITER_BUFFER_SIZE, WRITER_BUFFER_SIZE,
        tf::io::ZlibCompressionOptions::GZIP());
    TF_CHECK_OK(zwriter->Init());
    writer_.reset(zwriter);
  } else {
    writer_ = raw_file_;
  }
}

FastqWriter::~FastqWriter() {
  if (writer_) {
    TF_CHECK_OK(Close());
  }
}

tf::Status FastqWriter::Close() {
  if (!writer_)
    return tf::errors::FailedPrecondition(
        "Cannot close an already closed FastqWriter");
  // Close the file pointer we have been writing to.
  TF_RETURN_IF_ERROR(writer_->Close());
  writer_.reset();
  if (raw_file_) {
    // If this is a compressed file, the raw_file_ pointer is different and
    // also needs to be closed.
    if (isCompressed_) TF_RETURN_IF_ERROR(raw_file_->Close());
    raw_file_.reset();
  }
  return tf::Status::OK();
}

tf::Status FastqWriter::Write(
    const nucleus::genomics::v1::FastqRecord& record) {
  if (!writer_)
    return tf::errors::FailedPrecondition(
        "Cannot write to closed FASTQ stream.");
  string out = "@";
  tf::strings::StrAppend(&out, record.id());
  if (!record.description().empty()) {
    tf::strings::StrAppend(&out, " ", record.description());
  }
  tf::strings::StrAppend(&out, "\n", record.sequence(), "\n+\n",
                         record.quality(), "\n");
  TF_RETURN_IF_ERROR(writer_->Append(out));

  return tf::Status::OK();
}

}  // namespace nucleus
