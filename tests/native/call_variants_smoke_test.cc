// Phase 2 gate: smoke test for dv_tfrecord + dv_coreml.
//
// Tests TFRecord round-trip (write then read back) without needing a model.
// Core ML model loading is tested when the .mlpackage is available.

#include "gtest/gtest.h"
#include "deepvariant/native/tfrecord.h"

#include <filesystem>
#include <string>
#include <vector>

namespace deepvariant {
namespace {

// Write N records, read them back, verify content.
TEST(TFRecordRoundTrip, WriteAndReadBack) {
  auto tmp = std::filesystem::temp_directory_path() /
             "dv_tfrecord_test.tfrecord";
  const std::string path = tmp.string();

  std::vector<std::string> payloads = {
      "hello_world",
      std::string(10, '\x00'),         // all zeros
      std::string(1024, 'A'),          // 1 KB payload
      std::string("proto\x01\x02\x03"), // binary-looking data
  };

  // Write.
  {
    auto w = TFRecordWriter::New(path);
    ASSERT_NE(w, nullptr);
    for (const auto& p : payloads) {
      EXPECT_TRUE(w->WriteRecord(p));
    }
    EXPECT_TRUE(w->Flush());
    EXPECT_TRUE(w->Close());
  }

  // Read back.
  {
    auto r = TFRecordReader::New(path);
    ASSERT_NE(r, nullptr);
    for (size_t i = 0; i < payloads.size(); ++i) {
      ASSERT_TRUE(r->GetNext()) << "expected record " << i;
      EXPECT_EQ(r->record(), payloads[i]) << "record " << i << " mismatch";
    }
    EXPECT_FALSE(r->GetNext()) << "extra record found";
    r->Close();
  }

  std::filesystem::remove(path);
}

TEST(TFRecordRoundTrip, EmptyFile) {
  auto tmp = std::filesystem::temp_directory_path() /
             "dv_tfrecord_empty.tfrecord";
  const std::string path = tmp.string();
  {
    auto w = TFRecordWriter::New(path);
    ASSERT_NE(w, nullptr);
    w->Close();
  }
  {
    auto r = TFRecordReader::New(path);
    ASSERT_NE(r, nullptr);
    EXPECT_FALSE(r->GetNext());
  }
  std::filesystem::remove(path);
}

TEST(TFRecordRoundTrip, LargePayload) {
  auto tmp = std::filesystem::temp_directory_path() /
             "dv_tfrecord_large.tfrecord";
  const std::string path = tmp.string();
  // 100 * 221 * 7 * 4 bytes = one full pileup image as float32
  const size_t image_bytes = 100 * 221 * 7 * 4;
  std::string payload(image_bytes, '\x42');
  {
    auto w = TFRecordWriter::New(path);
    ASSERT_NE(w, nullptr);
    EXPECT_TRUE(w->WriteRecord(payload));
    w->Close();
  }
  {
    auto r = TFRecordReader::New(path);
    ASSERT_NE(r, nullptr);
    ASSERT_TRUE(r->GetNext());
    EXPECT_EQ(r->record().size(), image_bytes);
    EXPECT_EQ(r->record(), payload);
  }
  std::filesystem::remove(path);
}

}  // namespace
}  // namespace deepvariant
