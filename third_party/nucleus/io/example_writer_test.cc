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
#include "third_party/nucleus/io/example_writer.h"
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "tensorflow/core/lib/io/record_reader.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/test.h"

namespace nucleus {

TEST(ExampleWriterTest, TFErrorFilename) {
  auto output_record = std::string("test_output");
  EXPECT_DEATH(ExampleWriter("/tmp/out.csv"),
  "Unsupported file extension");
}

TEST(ExampleWriterTest, TFUpperCase) {
  auto output_record = std::string("test_output");
  ExampleWriter writer = ExampleWriter("/tmp/out.TFRECORD.GZ");
  EXPECT_TRUE(writer.Add(output_record));
}


TEST(ExampleWriterTest, NoParentDir) {
  // Tests for when no directory needs to be created.
  auto output_record = std::string("test_output");
  std::filesystem::current_path("/tmp");
  ExampleWriter writer = ExampleWriter("out.TFRECORD.GZ");
  EXPECT_TRUE(writer.Add(output_record));
}


TEST(ExampleWriterTest, TFWriterTest) {
  auto output_record = std::string("test_output");
  ExampleWriter writer = ExampleWriter("/tmp/out.tfrecord.gz");
  EXPECT_TRUE(writer.Add(output_record));
}

TEST(ExampleWriterTest, CompressionTypeForPathDetectsCodec) {
  EXPECT_EQ(CompressionTypeForPath("examples.tfrecord.snappy"), "SNAPPY");
  EXPECT_EQ(CompressionTypeForPath("examples.tfrecord.SNAPPY"), "SNAPPY");
  EXPECT_EQ(CompressionTypeForPath("examples-00000-of-00010.tfrecord.snappy"),
            "SNAPPY");
  EXPECT_EQ(CompressionTypeForPath("examples.tfrecord.gz"), "GZIP");
  EXPECT_EQ(CompressionTypeForPath("examples.tfrecord"), "GZIP");
}

TEST(ExampleWriterTest, SnappyRoundTrip) {
  // A ".snappy" path must produce a file a SNAPPY RecordReader can decode,
  // proving the writer's suffix-inferred codec agrees with the readers.
  const std::string path = "/tmp/dv_example_writer_snappy.tfrecord.snappy";
  const std::vector<std::string> records = {"alpha", "beta", "gamma"};
  {
    ExampleWriter writer(path);
    for (const std::string& r : records) EXPECT_TRUE(writer.Add(r));
    EXPECT_TRUE(writer.Close());
  }
  std::unique_ptr<tensorflow::RandomAccessFile> file;
  ASSERT_TRUE(
      tensorflow::Env::Default()->NewRandomAccessFile(path, &file).ok());
  tensorflow::io::RecordReaderOptions opts =
      tensorflow::io::RecordReaderOptions::CreateRecordReaderOptions("SNAPPY");
  tensorflow::io::RecordReader reader(file.get(), opts);
  tensorflow::uint64 offset = 0;
  tensorflow::tstring value;
  for (const std::string& expected : records) {
    ASSERT_TRUE(reader.ReadRecord(&offset, &value).ok());
    EXPECT_EQ(std::string(value.data(), value.size()), expected);
  }
}

TEST(ExampleWriterTest, GzipCompressionLevelIsApplied) {
  // A compressible payload, so level 0 (store) yields a strictly larger file
  // than level 9 -- proving the level reaches the zlib options.
  const std::string record(4096, 'A');
  auto write_at_level = [&record](const std::string& path, int level) {
    ExampleWriter writer(path, ExampleFormat::kAuto, level);
    for (int i = 0; i < 50; ++i) EXPECT_TRUE(writer.Add(record));
    EXPECT_TRUE(writer.Close());
  };
  const std::string path0 = "/tmp/dv_example_writer_l0.tfrecord.gz";
  const std::string path9 = "/tmp/dv_example_writer_l9.tfrecord.gz";
  write_at_level(path0, 0);
  write_at_level(path9, 9);
  tensorflow::uint64 size0 = 0, size9 = 0;
  ASSERT_TRUE(tensorflow::Env::Default()->GetFileSize(path0, &size0).ok());
  ASSERT_TRUE(tensorflow::Env::Default()->GetFileSize(path9, &size9).ok());
  EXPECT_GT(size0, size9);
}

}  // namespace nucleus
