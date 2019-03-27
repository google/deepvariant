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

#ifndef THIRD_PARTY_NUCLEUS_IO_TFRECORD_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_TFRECORD_WRITER_H_

#include <memory>
#include <string>

namespace tensorflow {
class WritableFile;
namespace io {
class RecordWriter;
}  // namespace io
}  // namespace tensorflow

namespace nucleus {

// A class for writing TFRecord files, designed for easy CLIF-wrapping
// for Python.  Loosely based on tensorflow/python/lib/io/py_record_writer.h
// An instance of this class is NOT safe for concurrent access by multiple
// threads.
class TFRecordWriter {
 public:
  // Create a TFRecordWriter.
  // Valid compression_types are "ZLIB", "GZIP", or "" (for none).
  // Returns nullptr on failure.
  static std::unique_ptr<TFRecordWriter> New(
      const std::string& filename, const std::string& compression_type);

  ~TFRecordWriter();

  // Returns true on success, false on error.
  bool WriteRecord(const std::string& record);

  // Returns true on success, false on error.
  bool Flush();

  // Close the file and release its resources.
  bool Close();

  // Disallow copy and assignment operations.
  TFRecordWriter(const TFRecordWriter& other) = delete;
  TFRecordWriter& operator=(const TFRecordWriter&) = delete;

 private:
  TFRecordWriter();

  // |writer_| has a non-owning pointer on |file_|, so destruct it first.
  std::unique_ptr<tensorflow::WritableFile> file_;
  std::unique_ptr<tensorflow::io::RecordWriter> writer_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_TFRECORD_WRITER_H_
