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

// Implementation of bed_writer.h
#include "third_party/nucleus/io/bed_writer.h"

#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/protos/bed.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/zlib_compression_options.h"
#include "third_party/nucleus/vendor/zlib_outputbuffer.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
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
// Writer for BED formats containing NGS reads.
//
// -----------------------------------------------------------------------------

StatusOr<std::unique_ptr<BedWriter>> BedWriter::ToFile(
    const string& bed_path, const nucleus::genomics::v1::BedHeader& header,
    const nucleus::genomics::v1::BedWriterOptions& options) {
  std::unique_ptr<tensorflow::WritableFile> fp;
  if (!tf::Env::Default()->NewWritableFile(bed_path.c_str(), &fp).ok()) {
    return tf::errors::Unknown(
        tf::strings::StrCat("Could not open bed_path ", bed_path));
  }
  bool isCompressed = EndsWith(bed_path, ".gz");
  auto writer = absl::WrapUnique(
      new BedWriter(std::move(fp), header, options, isCompressed));
  return std::move(writer);
}

BedWriter::BedWriter(std::unique_ptr<tensorflow::WritableFile> fp,
                     const nucleus::genomics::v1::BedHeader& header,
                     const nucleus::genomics::v1::BedWriterOptions& options,
                     const bool isCompressed)
    : header_(header),
      options_(options),
      raw_file_(std::move(fp)),
      isCompressed_(isCompressed) {
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

BedWriter::~BedWriter() {
  if (writer_) {
    TF_CHECK_OK(Close());
  }
}

tf::Status BedWriter::Close() {
  if (!writer_)
    return tf::errors::FailedPrecondition(
        "Cannot close an already closed BedWriter");
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

tf::Status BedWriter::Write(const nucleus::genomics::v1::BedRecord& record) {
  if (!writer_)
    return tf::errors::FailedPrecondition("Cannot write to closed BED stream.");
  int numFields = header_.num_fields();
  string out = "";
  absl::StrAppend(&out,
                  record.reference_name(),
                  "\t", record.start(),
                  "\t", record.end());
  if (numFields > 3) absl::StrAppend(&out, "\t", record.name());
  if (numFields > 4) absl::StrAppend(&out, "\t", record.score());
  if (numFields > 5) {
    switch (record.strand()) {
      case nucleus::genomics::v1::BedRecord::FORWARD_STRAND:
        absl::StrAppend(&out, "\t+");
        break;
      case nucleus::genomics::v1::BedRecord::REVERSE_STRAND:
        absl::StrAppend(&out, "\t-");
        break;
      case nucleus::genomics::v1::BedRecord::NO_STRAND:
        absl::StrAppend(&out, "\t.");
        break;
      default:
        return tf::errors::Unknown("Unknown strand encoding in the BED proto.");
    }
  }
  if (numFields > 7)
    absl::StrAppend(&out,
                    "\t", record.thick_start(),
                    "\t", record.thick_end());
  if (numFields > 8) absl::StrAppend(&out, "\t", record.item_rgb());
  if (numFields == 12)
    absl::StrAppend(&out,
                    "\t", record.block_count(),
                    "\t", record.block_sizes(),
                    "\t", record.block_starts());
  absl::StrAppend(&out, "\n");
  TF_RETURN_IF_ERROR(writer_->Append(out));

  return tf::Status::OK();
}

}  // namespace nucleus
