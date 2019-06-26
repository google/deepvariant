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

// Implementation of TextWriter class.
#include "third_party/nucleus/io/text_writer.h"

#include <stddef.h>
#include <string.h>
#include <sys/types.h>
#include <utility>

#include "absl/memory/memory.h"
#include "absl/strings/match.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "third_party/nucleus/io/hts_path.h"
#include "tensorflow/core/platform/logging.h"

namespace tf = tensorflow;

namespace {

// Write a string to an htslib file handle (compressed or not).
// Parallels hts_getline; oddly, no function like this is exposed by
// htslib.
tensorflow::Status hts_write(htsFile* hts_file, const char *str) {
  ssize_t str_len = strlen(str);
  ssize_t bytes_written;

  switch (hts_file->format.compression) {
    case no_compression:
      bytes_written = hwrite(hts_file->fp.hfile, str, str_len);
      break;
    case gzip:  // FALLTHROUGH_INTENDED
    case bgzf:
      bytes_written = bgzf_write(hts_file->fp.bgzf, str, str_len);
      break;
    default:
      return tf::errors::FailedPrecondition(
          "Unrecognized hts_file compression format");
  }

  if (bytes_written != str_len) {
    return tf::errors::DataLoss("Failure to write to htsFile.");
  }
  return tf::Status::OK();
}

}  // namespace



namespace nucleus {

StatusOr<std::unique_ptr<TextWriter>> TextWriter::ToFile(
    const string& path, CompressionPolicy compression) {
  const char* mode = compression == COMPRESS ? "wb" : "w";
  htsFile* fp = hts_open_x(path, mode);

  if (fp == nullptr) {
    return tf::errors::Unknown(
        "Could not open file for writing: ", path);
  } else {
    auto writer = absl::WrapUnique(new TextWriter(fp));
    return std::move(writer);
  }
}


StatusOr<std::unique_ptr<TextWriter>> TextWriter::ToFile(const string& path) {
  CompressionPolicy should_compress = (absl::EndsWith(path, ".gz") ?
                                       COMPRESS : NO_COMPRESS);
  return ToFile(path, should_compress);
}


TextWriter::TextWriter(htsFile* hts_file)
    : hts_file_(hts_file) {
  CHECK(hts_file_ != nullptr);
}

TextWriter::~TextWriter() {
  if (hts_file_) {
    TF_CHECK_OK(Close());
  }
}

tf::Status TextWriter::Write(const string& text) {
  if (hts_file_ == nullptr) {
    return tf::errors::FailedPrecondition(
        "Cannot write to a closed TextWriter");
  }
  return hts_write(hts_file_, text.c_str());
}


tf::Status TextWriter::Close() {
  if (!hts_file_) {
    return tf::errors::FailedPrecondition(
        "Cannot close an already closed file writer");
  }
  int hts_ok = hts_close(hts_file_);
  hts_file_ = nullptr;
  if (hts_ok < 0) {
    return tf::errors::Internal("hts_close() failed with return code ", hts_ok);
  }
  return tf::Status::OK();
}


}  // namespace nucleus
