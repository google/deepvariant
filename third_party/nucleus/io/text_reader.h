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

#ifndef THIRD_PARTY_NUCLEUS_IO_TEXT_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_TEXT_READER_H_

#include <memory>
#include <string>

#include "absl/memory/memory.h"
#include "htslib/hts.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "third_party/nucleus/core/status.h"

namespace nucleus {


// The TextReader class allows reading text from a (possibly compressed) file.
class TextReader {
 public:
  // Factory method to construct a TextReader.
  // File compression is determined from file magic (contents), not filename.
  static StatusOr<std::unique_ptr<TextReader>> FromFile(const string& path);

  // Destructor; closes the file, if it's still open.
  ~TextReader();

  // Reads a single line from the file.
  // Returns:
  //  - the string line (excluding trailing newline) if read is successful;
  //  - a status of ::nucleus::OutOfRange if at end-of-file;
  //  - otherwise, an appropriate error Status.
  StatusOr<string> ReadLine();

  // Explicitly closes the underlying file stream.
  ::nucleus::Status Close();

 private:
  // Private constructor.
  TextReader(htsFile* hts_file);

  // Underlying htslib file stream.
  htsFile* hts_file_;
};


}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_TEXT_READER_H_
