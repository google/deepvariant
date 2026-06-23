/*
 * Copyright 2024 Google LLC.
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
#ifndef THIRD_PARTY_NUCLEUS_IO_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_READER_H_

#include <cstdint>
#include <memory>
#include <string>
#include "absl/status/status.h"
#include "absl/strings/string_view.h"




namespace nucleus {

enum class ExampleFormat {
  kAuto = 0,  // Autodetect format by file extension.
  kTfRecord = 1,
  kBagz = 2,
};

// Returns the TfRecord compression codec implied by a file-name suffix:
// "SNAPPY" for a (case-insensitive) ".snappy" extension, otherwise "GZIP".
// Mirrors the reader-side detection in
// dv_utils.compression_type_for_examples_path so the file name is the single
// source of truth for both writing and reading.
std::string CompressionTypeForPath(absl::string_view path);

// Local writer for records, supports only a single file.
class ExampleWriter {
 public:
  // The TfRecord compression codec is inferred from the file-name suffix: a
  // ".snappy" extension selects SNAPPY, otherwise GZIP. compression_level only
  // applies to GZIP; a negative value uses the library default
  // (Z_DEFAULT_COMPRESSION).
  explicit ExampleWriter(absl::string_view path,
                         ExampleFormat format = ExampleFormat::kAuto,
                         int compression_level = -1);
  ~ExampleWriter();
  bool Add(absl::string_view value,
           absl::string_view chrom = {},
           int64_t pos = 0);
  bool Close();
  absl::Status status() { return status_; }

 private:
  class Impl;
  class TfRecordImpl;

  std::unique_ptr<Impl> impl_;
  absl::Status status_ = absl::OkStatus();
};

}  // namespace nucleus
#endif  // THIRD_PARTY_NUCLEUS_IO_READER_H_
