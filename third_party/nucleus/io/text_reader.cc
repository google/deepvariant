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

#include "third_party/nucleus/io/text_reader.h"

#include <stdlib.h>

#include <utility>

#include "absl/memory/memory.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/core/status.h"

namespace nucleus {

StatusOr<std::unique_ptr<TextReader>> TextReader::FromFile(const string& path) {
  htsFile* fp = hts_open_x(path, "r");

  if (fp == nullptr) {
    return ::nucleus::NotFound(
        absl::StrCat("Could not open ", path,
                     ". The file might not exist, or the format "
                     "detected by htslib might be incorrect."));
  } else {
    auto reader = absl::WrapUnique(new TextReader(fp));
    return std::move(reader);
  }
}

TextReader::~TextReader() {
  if (hts_file_) {
    NUCLEUS_CHECK_OK(Close());
  }
}

StatusOr<string> TextReader::ReadLine() {
  ::nucleus::Status status;
  string line;
  kstring_t k_line = {0, 0, nullptr};

  int ret = hts_getline(hts_file_, '\n', &k_line);
  if (ret == -1) {
    status = ::nucleus::OutOfRange("EOF");
  } else if (ret < 0) {
    status = ::nucleus::DataLoss("Failed to read text line");
  }

  if (k_line.s) {
    line = std::string(k_line.s);
    free(k_line.s);
  }
  if (status.ok()) {
    return line;
  } else {
    return status;
  }
}

::nucleus::Status TextReader::Close() {
  if (!hts_file_) {
    return ::nucleus::FailedPrecondition(
        "Cannot close an already closed file writer");
  }
  int hts_ok = hts_close(hts_file_);
  hts_file_ = nullptr;
  if (hts_ok < 0) {
    return ::nucleus::Internal(
        absl::StrCat("hts_close() failed with return code ", hts_ok));
  }
  return ::nucleus::Status();
}

TextReader::TextReader(htsFile* hts_file) : hts_file_(hts_file) {
  CHECK(hts_file_ != nullptr);
}

}  // namespace nucleus
