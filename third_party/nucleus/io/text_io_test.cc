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
#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/io/text_writer.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/core/statusor.h"

// -----------------------------------------------------------------------------
// Utility functions
namespace {

using nucleus::string;

const string FileContents(const string& path) {
  std::ifstream ifs(path, std::ifstream::binary);
  std::ostringstream oss;
  oss << ifs.rdbuf();
  return oss.str();
}

string MakeTempFileWithContents(const string& filename,
                                const string& contents) {
  string path = nucleus::MakeTempFile(filename);
  std::ofstream ofs(path, std::ofstream::binary);
  ofs << contents;
  return path;
}

}  // namespace


namespace nucleus {

// -----------------------------------------------------------------------------
// Test data

const char kHelloWorld[] = "Hello, world!";
const string kHelloWorldStr(kHelloWorld);

// Array generated via:
//  printf 'Hello, world!' | bgzip -c | hexdump -v -e '"0x" /1 "%02X" ", "'
const unsigned char kHelloWorldBGZF[] = {
  0x1F, 0x8B, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0x06, 0x00,
  0x42, 0x43, 0x02, 0x00, 0x28, 0x00, 0xF3, 0x48, 0xCD, 0xC9, 0xC9, 0xD7,
  0x51, 0x28, 0xCF, 0x2F, 0xCA, 0x49, 0x51, 0x04, 0x00, 0xE6, 0xC6, 0xE6,
  0xEB, 0x0D, 0x00, 0x00, 0x00, 0x1F, 0x8B, 0x08, 0x04, 0x00, 0x00, 0x00,
  0x00, 0x00, 0xFF, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1B, 0x00, 0x03,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
const string kHelloWorldBGZFStr(reinterpret_cast<const char*>(kHelloWorldBGZF),
                                sizeof(kHelloWorldBGZF));


// -----------------------------------------------------------------------------
// TextWriter tests

// Tests that we can write to an uncompressed file stream.
TEST(TextWriterTest, WritesUncompressedOutput) {
  string dest = MakeTempFile("uncompressed.txt");
  const auto writer = std::move(
      TextWriter::ToFile(dest, TextWriter::NO_COMPRESS).ValueOrDie());
  ::nucleus::Status status;
  status = writer->Write(kHelloWorld);
  EXPECT_EQ(::nucleus::Status(), status);
  status = writer->Close();
  EXPECT_EQ(::nucleus::Status(), status);

  EXPECT_EQ(kHelloWorldStr, FileContents(dest));
}

// Tests that we can write to a compresed file stream.
TEST(TextWriterTest, WritesCompressedOutput) {
  string dest = MakeTempFile("compressed.txt.gz");
  const auto writer = std::move(
      TextWriter::ToFile(dest, TextWriter::COMPRESS).ValueOrDie());
  ::nucleus::Status status;
  status = writer->Write(kHelloWorld);
  EXPECT_EQ(::nucleus::Status(), status);
  status = writer->Close();
  EXPECT_EQ(::nucleus::Status(), status);

  EXPECT_EQ(kHelloWorldBGZFStr, FileContents(dest));
}

// Tests that we get compression when filename is *.gz.
TEST(TextWriterTest, UsesCompressionWhenExtensionIsGz) {
  string destGz = MakeTempFile("test_file.txt.gz");
  const auto writer = std::move(TextWriter::ToFile(destGz).ValueOrDie());
  ::nucleus::Status status;
  status = writer->Write(kHelloWorld);
  EXPECT_EQ(::nucleus::Status(), status);
  status = writer->Close();
  EXPECT_EQ(::nucleus::Status(), status);

  EXPECT_EQ(kHelloWorldBGZFStr, FileContents(destGz));
}

// Tests that we get  no  compression when filename is NOT *.gz.
TEST(TextWriterTest, NoCompressionWhenExtensionIsNotGz) {
  string dest = MakeTempFile("test_file.txt");
  const auto writer = std::move(TextWriter::ToFile(dest).ValueOrDie());
  ::nucleus::Status status;
  status = writer->Write(kHelloWorld);
  EXPECT_EQ(::nucleus::Status(), status);
  status = writer->Close();
  EXPECT_EQ(::nucleus::Status(), status);

  EXPECT_EQ(kHelloWorld, FileContents(dest));
}

// -----------------------------------------------------------------------------
// TextReader tests

// Tests that we can read from an uncompressed file stream.
TEST(TextReaderTest, ReadsUncompressedFile) {
  string path =
      MakeTempFileWithContents("uncompressed-for-reader.txt", kHelloWorldStr);
  const auto reader = std::move(TextReader::FromFile(path).ValueOrDie());

  StatusOr<string> rv;
  rv = reader->ReadLine();
  EXPECT_TRUE(rv.ok());
  EXPECT_EQ(rv.ValueOrDie(), kHelloWorldStr);
  rv = reader->ReadLine();
  EXPECT_TRUE(::nucleus::IsOutOfRange(rv.status()));
}


// Tests that we can read from an uncompressed file stream.
TEST(TextReaderTest, ReadsCompressedFile) {
  string path =
      MakeTempFileWithContents("compressed-for-reader.bin", kHelloWorldBGZFStr);
  const auto reader = std::move(TextReader::FromFile(path).ValueOrDie());

  StatusOr<string> rv;
  rv = reader->ReadLine();
  EXPECT_TRUE(rv.ok());
  EXPECT_EQ(rv.ValueOrDie(), kHelloWorldStr);
  rv = reader->ReadLine();
  EXPECT_TRUE(::nucleus::IsOutOfRange(rv.status()));
}


}  // namespace nucleus
