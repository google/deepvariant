/*
 * Copyright 2018 Google Inc.
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

// Testing C++ utilities.
#ifndef THIRD_PARTY_NUCLEUS_TESTING_TEST_UTILS_H_
#define THIRD_PARTY_NUCLEUS_TESTING_TEST_UTILS_H_

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>


#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/io/record_reader.h"
#include "tensorflow/core/lib/io/record_writer.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

using tensorflow::StringPiece;
using tensorflow::string;
using tensorflow::uint64;

constexpr char kBioTFCoreTestDataDir[] = "third_party/nucleus/testdata";

// N.B. this will be set to "" in OSS.
constexpr char kDefaultWorkspace[] = "";


// Simple getter for test files in the right testdata path.
// This uses JoinPath, so no leading or trailing "/" are necessary.
string GetTestData(
    tensorflow::StringPiece path,
    tensorflow::StringPiece test_data_dir = kBioTFCoreTestDataDir);

// Returns a path to a temporary file with filename in the appropriate test
// directory.
string MakeTempFile(tensorflow::StringPiece filename);

// Reads all of the records from path into a vector of parsed Proto. Path
// must point to a TFRecord formatted file.
template <typename Proto>
std::vector<Proto> ReadProtosFromTFRecord(tensorflow::StringPiece path) {
  tensorflow::Env* env = tensorflow::Env::Default();
  std::unique_ptr<tensorflow::RandomAccessFile> read_file;
  TF_CHECK_OK(env->NewRandomAccessFile(path.ToString(), &read_file));
  tensorflow::io::RecordReader reader(read_file.get());
  std::vector<Proto> results;
  uint64 offset = 0;
  string data;
  while (reader.ReadRecord(&offset, &data).ok()) {
    Proto proto;
    CHECK(proto.ParseFromString(data)) << "Failed to parse proto";
    results.push_back(proto);
  }
  return results;
}

// Writes all `protos` to a TFRecord formatted file.
template <typename Proto>
void WriteProtosToTFRecord(const std::vector<Proto>& protos,
                           tensorflow::StringPiece output_path) {
  std::unique_ptr<tensorflow::WritableFile> file;
  TF_CHECK_OK(tensorflow::Env::Default()->NewWritableFile(
      output_path.ToString(), &file));
  tensorflow::io::RecordWriter record_writer(file.get());
  for (const Proto& proto : protos) {
    TF_CHECK_OK(record_writer.WriteRecord(proto.SerializeAsString()));
  }
}

// Creates a vector of ContigInfos with specified `names` and `positions`
// representing pos_in_fast. `names` and `positions` need to have the same
// number of elements.
std::vector<nucleus::genomics::v1::ContigInfo> CreateContigInfos(
    const std::vector<string>& names, const std::vector<int>& positions);

// A matcher to help us do Pointwise double expectations with a provided
// precision via DoubleNear.
MATCHER_P(PointwiseDoubleNear, abs_error, "") {
  using RhsType = decltype(::testing::get<1>(arg));
  ::testing::Matcher<RhsType> matcher =
      ::testing::DoubleNear(::testing::get<1>(arg), abs_error);
  return matcher.Matches(::testing::get<0>(arg));
}

// A matcher to test if a floating point value is finite.
MATCHER(IsFinite, "") { return std::isfinite(arg); }

// Adapter to extract an iterable into a vector for examination in test code
// from a StatusOr<std::shared_ptr<Iterable<Record>>>.
template <class Record>
std::vector<Record> as_vector(
    const StatusOr<std::shared_ptr<Iterable<Record>>>& it) {
  TF_CHECK_OK(it.status());
  return as_vector(it.ValueOrDie());
}

// Adapter to extract an iterable into a vector for examination in test code.
template<class Record>
std::vector<Record> as_vector(const std::shared_ptr<Iterable<Record>>& it) {
  std::vector<Record> records;
  for (const StatusOr<Record*> value_status : it) {
    records.push_back(*value_status.ValueOrDie());
  }
  return records;
}

// Creates a test Read.
//
// The read has reference_name chr, start of start, aligned_sequence of bases,
// and cigar element parsed from cigar_elements, which is vector of standard
// CIGAR element string values like {"5M", "2I", "3M"} which is 5 bp matches,
// 2 bp insertion, and 3 bp matches. The read has base qualities set to 30 and
// a mapping quality of 90.
::nucleus::genomics::v1::Read MakeRead(
    const string& chr, int start, const string& bases,
    const std::vector<string>& cigar_elements);

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_TESTING_TEST_UTILS_H_
