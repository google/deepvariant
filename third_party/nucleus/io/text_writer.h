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

#ifndef THIRD_PARTY_NUCLEUS_IO_TEXT_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_TEXT_WRITER_H_

#include <memory>
#include <string>

#include "htslib/hts.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

// TextWriter is a class allowing writing text to a (possibly compressed) file
// stream
class TextWriter {
 public:  // Types.
  enum CompressionPolicy {
    NO_COMPRESS = false,
    COMPRESS = true,
  };

 public:
  // Factory function allowing explicit choice of whether to use compression.
  static StatusOr<std::unique_ptr<TextWriter>> ToFile(
      const string& path, CompressionPolicy compression);

  // Factory function that uses compression if the filename ends in ".gz".
  static StatusOr<std::unique_ptr<TextWriter>> ToFile(const string& path);

  // Destructor; closes the stream if still open.
  ~TextWriter();

  // Write a string to the file stream.
  ::nucleus::Status Write(const string& text);

  // Close the underlying file stream.
  ::nucleus::Status Close();

 private:
  // Private constructor.
  TextWriter(htsFile* hts_file);

  // Underlying htslib file stream.
  htsFile* hts_file_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_TEXT_WRITER_H_
