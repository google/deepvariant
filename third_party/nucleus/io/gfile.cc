/*
 * Copyright 2019 Google LLC.
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

#include "third_party/nucleus/io/gfile.h"

#include "absl/memory/memory.h"
#include "tensorflow/core/lib/io/buffered_inputstream.h"
#include "tensorflow/core/lib/io/random_inputstream.h"
#include "tensorflow/core/platform/env.h"

namespace nucleus {

bool Exists(const std::string& filename) {
  // FileExists sets s to tensorflow::error::NOT_FOUND if it doesn't exist.
  tensorflow::Status s = tensorflow::Env::Default()->FileExists(filename);
  return s.ok();
}

std::vector<std::string> Glob(const std::string& pattern) {
  std::vector<std::string> results;
  tensorflow::Status s = tensorflow::Env::Default()->GetMatchingPaths(
      pattern, &results);
  return results;
}

ReadableFile::ReadableFile() {}

std::unique_ptr<ReadableFile> ReadableFile::New(const std::string& filename) {
  std::unique_ptr<tensorflow::RandomAccessFile> file;
  tensorflow::Status status =
      tensorflow::Env::Default()->NewRandomAccessFile(filename, &file);
  if (!status.ok()) {
    return nullptr;
  }

  size_t buffer_size = 512 * 1024;

  std::unique_ptr<tensorflow::io::RandomAccessInputStream> input_stream(
      new tensorflow::io::RandomAccessInputStream(
          file.release(), true /* owns_file */));
  std::unique_ptr<tensorflow::io::BufferedInputStream> buffered_input_stream(
      new tensorflow::io::BufferedInputStream(
          input_stream.release(), buffer_size, true /* owns_input_stream */));

  auto f = absl::WrapUnique<ReadableFile>(new ReadableFile);
  f->stream_ = std::move(buffered_input_stream);

  return f;
}

ReadableFile::~ReadableFile() {
}

bool ReadableFile::Readline(std::string* s) {
  *s = stream_->ReadLineAsString();
  return s->length() > 0;
}

void ReadableFile::Close() {
  stream_ = nullptr;
}

WritableFile::WritableFile() {}

std::unique_ptr<WritableFile> WritableFile::New(const std::string& filename) {
  std::unique_ptr<tensorflow::WritableFile> file;

  tensorflow::Status s = tensorflow::Env::Default()->NewWritableFile(
      filename, &file);

  if (!s.ok()) {
    return nullptr;
  }

  auto f = absl::WrapUnique<WritableFile>(new WritableFile);
  f->file_ = std::move(file);

  return f;
}

bool WritableFile::Write(const std::string& s) {
  tensorflow::Status status = file_->Append(s);
  return status.ok();
}

void WritableFile::Close() {
  file_ = nullptr;
}

WritableFile::~WritableFile() {
}

}  // namespace nucleus

