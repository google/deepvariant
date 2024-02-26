/*
 * Copyright 2022 Google LLC.
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

#include "third_party/nucleus/io/variant_reader.h"

#include <string>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/io/tfrecord_writer.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"

namespace {

using nucleus::EqualsProto;
using testing::Pointee;

nucleus::genomics::v1::Variant VariantProto(const std::string& ref_name,
                                            int start, int end = 0) {
  nucleus::genomics::v1::Variant v;
  v.set_reference_name(ref_name);
  v.set_start(start);
  if (end != 0) {
    v.set_end(end);
  }

  return v;
}

std::string VariantStr(const std::string& ref_name, int start, int end = 0) {
  return VariantProto(ref_name, start, end).SerializeAsString();
}

TEST(IndexedReaderTest, EmptyShards) {
  std::string path_a = absl::StrCat(::testing::TempDir(), "/", "a.gz");
  std::string path_b = absl::StrCat(::testing::TempDir(), "/", "b.gz");
  auto writer_a = nucleus::TFRecordWriter::New(path_a, "GZIP");
  auto writer_b = nucleus::TFRecordWriter::New(path_b, "GZIP");
  writer_a->Close();
  writer_b->Close();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map = {{}};

  auto reader =
      nucleus::ShardedVariantReader::Open({path_a, path_b}, contig_index_map);
  EXPECT_EQ(reader->GetAndReadNext().variant, nullptr);
}

TEST(IndexedReaderTest, ReadsRecordsSingleShard) {
  std::string path_a = absl::StrCat(::testing::TempDir(), "/", "a.gz");
  auto writer_a = nucleus::TFRecordWriter::New(path_a, "GZIP");
  writer_a->WriteRecord(VariantStr("ref_a", 1));
  writer_a->WriteRecord(VariantStr("ref_b", 1));
  writer_a->WriteRecord(VariantStr("ref_c", 1));
  writer_a->Close();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map = {
      {"ref_a", 0}, {"ref_b", 1}, {"ref_c", 2}};
  auto reader = nucleus::ShardedVariantReader::Open({path_a}, contig_index_map);
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_b", 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_c", 1))));
}

TEST(IndexedReaderTest, RecordsAreGloballyOrdered) {
  std::string path_a = absl::StrCat(::testing::TempDir(), "/", "a.gz");
  auto writer_a = nucleus::TFRecordWriter::New(path_a, "GZIP");
  writer_a->WriteRecord(VariantStr("ref_a", 2));
  writer_a->WriteRecord(VariantStr("ref_c", 100));
  writer_a->WriteRecord(VariantStr("ref_d", 1));
  writer_a->Close();

  std::string path_b = absl::StrCat(::testing::TempDir(), "/", "b.gz");
  auto writer_b = nucleus::TFRecordWriter::New(path_b, "GZIP");
  writer_b->WriteRecord(VariantStr("ref_a", 1));
  writer_b->WriteRecord(VariantStr("ref_b", 1));
  writer_b->WriteRecord(VariantStr("ref_d", 2));
  writer_b->WriteRecord(VariantStr("ref_d", 3));
  writer_b->WriteRecord(VariantStr("ref_d", 4));
  writer_b->Close();

  std::string path_empty = absl::StrCat(::testing::TempDir(), "/", "empty.gz");
  auto writer_c = nucleus::TFRecordWriter::New(path_empty, "GZIP");
  writer_c->Close();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map = {
      {"ref_a", 0}, {"ref_b", 1}, {"ref_c", 2}, {"ref_d", 3}};
  auto reader = nucleus::ShardedVariantReader::Open(
      {path_a, path_b, path_empty}, contig_index_map);

  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 2))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_b", 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_c", 100))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_d", 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_d", 2))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_d", 3))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_d", 4))));
}

TEST(IndexedReaderTest, RecordsAreGloballyOrderedWithinContigs) {
  std::string path_a = absl::StrCat(::testing::TempDir(), "/", "a.gz");
  auto writer_a = nucleus::TFRecordWriter::New(path_a, "GZIP");
  writer_a->WriteRecord(VariantStr("ref_a", 2, 3));
  writer_a->WriteRecord(VariantStr("ref_a", 2, 5));
  writer_a->WriteRecord(VariantStr("ref_a", 3, 4));
  writer_a->WriteRecord(VariantStr("ref_b", 1, 1));
  writer_a->Close();

  std::string path_b = absl::StrCat(::testing::TempDir(), "/", "b.gz");
  auto writer_b = nucleus::TFRecordWriter::New(path_b, "GZIP");
  writer_b->WriteRecord(VariantStr("ref_a", 2, 2));
  writer_b->WriteRecord(VariantStr("ref_a", 2, 6));
  writer_b->WriteRecord(VariantStr("ref_a", 3, 1));
  writer_b->WriteRecord(VariantStr("ref_b", 2, 1));
  writer_b->Close();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map = {{"ref_a", 0},
                                                                 {"ref_b", 1}};
  auto reader =
      nucleus::ShardedVariantReader::Open({path_a, path_b}, contig_index_map);

  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 2, 2))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 2, 3))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 2, 5))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 2, 6))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 3, 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_a", 3, 4))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_b", 1, 1))));
  EXPECT_THAT(reader->GetAndReadNext().variant,
              Pointee(EqualsProto(VariantProto("ref_b", 2, 1))));
}

TEST(IndexedReaderTest, ReturnsContigIndex) {
  std::string path_a = absl::StrCat(::testing::TempDir(), "/", "a");
  auto writer_a = nucleus::TFRecordWriter::New(path_a, "");
  writer_a->WriteRecord(VariantStr("ref_a", 1));
  writer_a->WriteRecord(VariantStr("ref_b", 2));
  writer_a->Close();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map = {{"ref_a", 1},
                                                                 {"ref_b", 2}};
  auto reader = nucleus::ShardedVariantReader::Open({path_a}, contig_index_map);
  EXPECT_THAT(reader->GetAndReadNext().contig_map_index, 1);
  EXPECT_THAT(reader->GetAndReadNext().contig_map_index, 2);
}

}  // namespace
