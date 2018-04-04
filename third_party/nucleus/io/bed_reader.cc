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
 */

// Implementation of bed_reader.h
#include "third_party/nucleus/io/bed_reader.h"

#include "absl/strings/str_split.h"
#include "third_party/nucleus/protos/bed.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/zlib_compression_options.h"
#include "third_party/nucleus/vendor/zlib_inputstream.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/io/buffered_inputstream.h"
#include "tensorflow/core/lib/io/compression.h"
#include "tensorflow/core/lib/io/random_inputstream.h"
#include "tensorflow/core/lib/strings/numbers.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/file_system.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

namespace tf = tensorflow;

using tensorflow::int32;
using tensorflow::int64;
using tensorflow::string;

// 256 KB read buffer.
constexpr int READER_BUFFER_SIZE = 256 * 1024;

// BED-specific attributes.
constexpr char BED_COMMENT_CHAR = '#';

// -----------------------------------------------------------------------------
//
// Reader for BED format data.
//
// -----------------------------------------------------------------------------

namespace {

bool ValidNumBedFields(const int fields) {
  return (fields == 3 || fields == 4 || fields == 5 || fields == 6 ||
          fields == 8 || fields == 9 || fields == 12);
}

// Read the next non-comment line.
tf::Status NextNonCommentLine(
    const std::unique_ptr<tf::io::BufferedInputStream>& instream,
    string* line) {
  TF_RETURN_IF_ERROR(instream->ReadLine(line));
  while ((*line)[0] == BED_COMMENT_CHAR) {
    TF_RETURN_IF_ERROR(instream->ReadLine(line));
  }
  return tf::Status::OK();
}

tf::Status ConvertToPb(const string& line, const int desiredNumFields,
                       int* numTokensSeen,
                       nucleus::genomics::v1::BedRecord* record) {
  CHECK(record != nullptr) << "BED record cannot be null";
  record->Clear();

  std::vector<string> tokens = absl::StrSplit(line, '\t');
  int numTokens = static_cast<int>(tokens.size());
  *numTokensSeen = numTokens;
  if (!ValidNumBedFields(numTokens)) {
    return tf::errors::Unknown("BED record has invalid number of fields");
  }
  int numFields =
      desiredNumFields == 0 ? numTokens : std::min(numTokens, desiredNumFields);
  int64 int64Value;
  record->set_reference_name(tokens[0]);
  tf::strings::safe_strto64(tokens[1], &int64Value);
  record->set_start(int64Value);
  tf::strings::safe_strto64(tokens[2], &int64Value);
  record->set_end(int64Value);
  if (numFields > 3) record->set_name(tokens[3]);
  if (numFields > 4) {
    double value;
    tf::strings::safe_strtod(tokens[4].c_str(), &value);
    record->set_score(value);
  }
  if (numFields > 5) {
    if (tokens[5] == "+")
      record->set_strand(nucleus::genomics::v1::BedRecord::FORWARD_STRAND);
    else if (tokens[5] == "-")
      record->set_strand(nucleus::genomics::v1::BedRecord::REVERSE_STRAND);
    else if (tokens[5] == ".")
      record->set_strand(nucleus::genomics::v1::BedRecord::NO_STRAND);
    else
      return tf::errors::Unknown("Invalid BED record with unknown strand");
  }
  if (numFields > 7) {
    tf::strings::safe_strto64(tokens[6], &int64Value);
    record->set_thick_start(int64Value);
    tf::strings::safe_strto64(tokens[7], &int64Value);
    record->set_thick_end(int64Value);
  }
  if (numFields > 8) record->set_item_rgb(tokens[8]);
  if (numFields >= 12) {
    int32 int32Value;
    tf::strings::safe_strto32(tokens[9], &int32Value);
    record->set_block_count(int32Value);
    record->set_block_sizes(tokens[10]);
    record->set_block_starts(tokens[11]);
  }

  return tf::Status::OK();
}

// Peeks into the path to the first BED record and returns the number of fields
// in the record.
// NOTE: This is quite heavyweight. Reading upon initialization and then
// rewinding the stream to 0 is a nicer solution, but currently has a memory
// leak in the compressed stream reset implementation.
tf::Status GetNumFields(const string& path, bool isCompressed, int* numFields) {
  std::unique_ptr<tensorflow::RandomAccessFile> sp;
  tf::Status status =
      tf::Env::Default()->NewRandomAccessFile(path.c_str(), &sp);
  if (!status.ok()) {
    return tf::errors::NotFound(tf::strings::StrCat("Could not open ", path));
  }
  tensorflow::RandomAccessFile* fp = sp.release();
  std::unique_ptr<tensorflow::io::BufferedInputStream> bi;
  string line;
  if (isCompressed) {
    std::unique_ptr<tensorflow::io::RandomAccessInputStream> fs;
    std::unique_ptr<tensorflow::io::ZlibInputStream> zs;
    fs.reset(new tf::io::RandomAccessInputStream(fp));
    zs.reset(new tf::io::ZlibInputStream(
        fs.get(), READER_BUFFER_SIZE, READER_BUFFER_SIZE,
        tf::io::ZlibCompressionOptions::GZIP()));
    bi.reset(new tf::io::BufferedInputStream(zs.get(), READER_BUFFER_SIZE));
    TF_RETURN_IF_ERROR(NextNonCommentLine(bi, &line));
    bi.reset();
    zs.reset();
    fs.reset();
  } else {
    bi.reset(new tf::io::BufferedInputStream(fp, READER_BUFFER_SIZE));
    TF_RETURN_IF_ERROR(NextNonCommentLine(bi, &line));
    bi.reset();
  }
  delete fp;
  std::vector<string> tokens = absl::StrSplit(line, '\t');
  *numFields = static_cast<int>(tokens.size());
  return tf::Status::OK();
}
}  // namespace

// Iterable class for traversing all BED records in the file.
class BedFullFileIterable : public BedIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::BedRecord* out) override;

  // Constructor is invoked via BedReader::Iterate.
  BedFullFileIterable(const BedReader* reader);
  ~BedFullFileIterable() override;
};

StatusOr<std::unique_ptr<BedReader>> BedReader::FromFile(
    const string& bed_path,
    const nucleus::genomics::v1::BedReaderOptions& options) {
  int numFieldsInBed;
  TF_RETURN_IF_ERROR(
      GetNumFields(bed_path,
                   options.compression_type() ==
                       nucleus::genomics::v1::BedReaderOptions::GZIP,
                   &numFieldsInBed));
  nucleus::genomics::v1::BedHeader header;
  header.set_num_fields(numFieldsInBed);
  // Ensure options are valid.
  if (options.num_fields() != 0 && (options.num_fields() > numFieldsInBed ||
                                    !ValidNumBedFields(options.num_fields()))) {
    return tf::errors::InvalidArgument(
        "Invalid requested number of fields to parse");
  }
  std::unique_ptr<tensorflow::RandomAccessFile> fp;
  tf::Status status =
      tf::Env::Default()->NewRandomAccessFile(bed_path.c_str(), &fp);
  if (!status.ok()) {
    return tf::errors::NotFound(
        tf::strings::StrCat("Could not open ", bed_path));
  }
  return std::unique_ptr<BedReader>(
      new BedReader(fp.release(), options, header));
}

BedReader::BedReader(tensorflow::RandomAccessFile* fp,
                     const nucleus::genomics::v1::BedReaderOptions& options,
                     const nucleus::genomics::v1::BedHeader& header)
    : options_(options), header_(header), src_(fp) {
  if (options.compression_type() ==
      nucleus::genomics::v1::BedReaderOptions::GZIP) {
    file_stream_.reset(new tf::io::RandomAccessInputStream(src_));
    zlib_stream_.reset(new tf::io::ZlibInputStream(
        file_stream_.get(), READER_BUFFER_SIZE, READER_BUFFER_SIZE,
        tf::io::ZlibCompressionOptions::GZIP()));
    buffered_inputstream_.reset(new tf::io::BufferedInputStream(
        zlib_stream_.get(), READER_BUFFER_SIZE));
  } else {
    buffered_inputstream_.reset(
        new tf::io::BufferedInputStream(src_, READER_BUFFER_SIZE));
  }
}

BedReader::~BedReader() {
  if (src_) {
    TF_CHECK_OK(Close());
  }
}

tf::Status BedReader::Close() {
  if (src_ == nullptr) {
    return tf::errors::FailedPrecondition("BedReader already closed");
  }
  buffered_inputstream_.reset();
  zlib_stream_.reset();
  file_stream_.reset();
  delete src_;
  src_ = nullptr;
  return tf::Status::OK();
}

// Ensures the number of fields is consistent across all records in the BED.
tf::Status BedReader::Validate(const int numTokens) const {
  if (header_.num_fields() != numTokens) {
    return tf::errors::Unknown(
        "Invalid BED with varying number of fields in file");
  }
  return tf::Status::OK();
}

StatusOr<std::shared_ptr<BedIterable>> BedReader::Iterate() const {
  if (src_ == nullptr)
    return tf::errors::FailedPrecondition("Cannot Iterate a closed BedReader.");
  return StatusOr<std::shared_ptr<BedIterable>>(
      MakeIterable<BedFullFileIterable>(this));
}

// Iterable class definitions.
StatusOr<bool> BedFullFileIterable::Next(
    nucleus::genomics::v1::BedRecord* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  const BedReader* bed_reader = static_cast<const BedReader*>(reader_);
  string line;
  tf::Status status = NextNonCommentLine(bed_reader->Stream(), &line);
  if (tf::errors::IsOutOfRange(status)) {
    return false;
  } else if (!status.ok()) {
    return status;
  }
  int numTokens;
  TF_RETURN_IF_ERROR(
      ConvertToPb(line, bed_reader->Options().num_fields(), &numTokens, out));
  TF_RETURN_IF_ERROR(bed_reader->Validate(numTokens));
  return true;
}

BedFullFileIterable::~BedFullFileIterable() {}

BedFullFileIterable::BedFullFileIterable(const BedReader* reader)
    : Iterable(reader) {}

}  // namespace nucleus
