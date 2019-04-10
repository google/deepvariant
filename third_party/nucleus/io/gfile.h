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

#ifndef THIRD_PARTY_NUCLEUS_IO_GFILE_H_
#define THIRD_PARTY_NUCLEUS_IO_GFILE_H_

#include <memory>
#include <string>
#include <vector>

namespace tensorflow {
class WritableFile;
namespace io {
class BufferedInputStream;
}  // namespace io
}  // namespace tensorflow

namespace nucleus {

// Return whether or not filename exists as a file.
bool Exists(const std::string& filename);

// Return all files matching the shell-style file glob.
std::vector<std::string> Glob(const std::string& pattern);

class ReadableFile {
 public:
  static std::unique_ptr<ReadableFile> New(const std::string& filename);
  ~ReadableFile();

  // Reads the next line into *s, and returns true if that went ok.
  bool Readline(std::string* s);

  void Close();

  // This no-op function is needed only for Python context manager support.
  void PythonEnter() const {}

 private:
  ReadableFile();

  std::unique_ptr<tensorflow::io::BufferedInputStream> stream_;
};

class WritableFile {
 public:
  static std::unique_ptr<WritableFile> New(const std::string& filename);
  ~WritableFile();

  bool Write(const std::string& s);
  void Close();

  // This no-op function is needed only for Python context manager support.
  void PythonEnter() const {}

 private:
  WritableFile();

  std::unique_ptr<tensorflow::WritableFile> file_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_GFILE_H_
