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

#include "third_party/nucleus/io/tfrecord_reader.h"

#include <memory>

#include "absl/log/log.h"
#include "absl/memory/memory.h"
#include "tensorflow/core/lib/io/record_reader.h"

namespace nucleus {

TFRecordReader::TFRecordReader() {}

std::unique_ptr<TFRecordReader> TFRecordReader::New(
    const std::string& filename, const std::string& compression_type) {
  std::unique_ptr<tensorflow::RandomAccessFile> file;
  tensorflow::Status s =
      tensorflow::Env::Default()->NewRandomAccessFile(filename, &file);
  if (!s.ok()) {
    LOG(ERROR) << s;
    return nullptr;
  }

  auto reader = absl::WrapUnique<TFRecordReader>(new TFRecordReader);
  reader->offset_ = 0;
  reader->file_ = std::move(file);

  tensorflow::io::RecordReaderOptions options =
      tensorflow::io::RecordReaderOptions::CreateRecordReaderOptions(
          compression_type);
  options.buffer_size = 16 * 1024 * 1024;
  reader->reader_ = std::make_unique<tensorflow::io::RecordReader>(
      reader->file_.get(), options);

  return reader;
}

TFRecordReader::~TFRecordReader() {
}

bool TFRecordReader::GetNext() {
  if (reader_ == nullptr) {
    return false;
  }

  tensorflow::Status s = reader_->ReadRecord(&offset_, &record_);

  return s.ok();
}

void TFRecordReader::Close() {
  reader_ = nullptr;
  file_ = nullptr;
}

}  // namespace nucleus
