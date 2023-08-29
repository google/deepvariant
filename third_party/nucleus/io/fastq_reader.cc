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
 */

// Implementation of fastq_reader.h
#include "third_party/nucleus/io/fastq_reader.h"

#include <stddef.h>

#include <utility>

#include "absl/strings/string_view.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/fastq.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

using absl::string_view;
using nucleus::genomics::v1::FastqReaderOptions;
using nucleus::genomics::v1::FastqRecord;

// For validation of the FASTQ format.
constexpr char HEADER_SYMBOL = '@';
constexpr char SEQUENCE_AND_QUALITY_SEPARATOR_SYMBOL = '+';

// -----------------------------------------------------------------------------
//
// Reader for FASTQ formats containing NGS reads.
//
// -----------------------------------------------------------------------------

namespace {

// TODO: get rid of pessimizing string_view -> string conversions
// once our OSS dependencies are updated.
::nucleus::Status ConvertToPb(const string_view header,
                              const string_view sequence, const string_view pad,
                              const string_view quality,
                              nucleus::genomics::v1::FastqRecord* record) {
  CHECK(record != nullptr) << "FASTQ record cannot be null";
  if (header.empty() || header[0] != HEADER_SYMBOL || pad.empty() ||
      pad[0] != SEQUENCE_AND_QUALITY_SEPARATOR_SYMBOL || sequence.empty() ||
      sequence.length() != quality.length()) {
    return ::nucleus::DataLoss("Invalid FASTQ record");
  }
  record->Clear();
  size_t spaceix = header.find(' ');
  if (spaceix == string::npos) {
    // No space found; ID is full string after delimiter.
    record->set_id(string(header.substr(1)));
  } else {
    // ID is the string from delimiter up to the first space.
    record->set_id(string(header.substr(1, spaceix - 1)));
    // Description is the string after the first space.
    record->set_description(string(header.substr(spaceix + 1)));
  }
  record->set_sequence(string(sequence));
  record->set_quality(string(quality));
  return ::nucleus::Status();
}
}  // namespace

// Iterable class for traversing all FASTQ records in the file.
class FastqFullFileIterable : public FastqIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::FastqRecord* out) override;

  // Constructor is invoked via FastqReader::Iterate.
  FastqFullFileIterable(const FastqReader* reader);
  ~FastqFullFileIterable() override;
};

StatusOr<std::unique_ptr<FastqReader>> FastqReader::FromFile(
    const string& fastq_path,
    const nucleus::genomics::v1::FastqReaderOptions& options) {
  StatusOr<std::unique_ptr<TextReader>> textreader_or =
      TextReader::FromFile(fastq_path);
  NUCLEUS_RETURN_IF_ERROR(textreader_or.status());
  return std::unique_ptr<FastqReader>(
      new FastqReader(std::move(textreader_or.ValueOrDie()), options));
}

FastqReader::FastqReader(std::unique_ptr<TextReader> text_reader,
                         const FastqReaderOptions& options)
    : options_(options), text_reader_(std::move(text_reader)) {}

FastqReader::~FastqReader() {
  if (text_reader_) {
    NUCLEUS_CHECK_OK(Close());
  }
}

::nucleus::Status FastqReader::Close() {
  if (!text_reader_) {
    return ::nucleus::FailedPrecondition("FastqReader already closed");
  }
  // Close the file pointer.
  ::nucleus::Status close_status = text_reader_->Close();
  text_reader_ = nullptr;
  return close_status;
}

::nucleus::Status FastqReader::Next(string* header, string* sequence,
                                    string* pad, string* quality) const {
  // Read the four lines, returning early if we are at the end of the stream or
  // the record is truncated.
  StatusOr<string> header_or, sequence_or, pad_or, quality_or;

  header_or = text_reader_->ReadLine();
  if (!header_or.ok()) {
    if (::nucleus::IsOutOfRange(header_or.status())) {
      return header_or.status();
    } else {
      goto data_loss;
    }
  }
  sequence_or = text_reader_->ReadLine();
  if (!sequence_or.ok()) {
    goto data_loss;
  }

  pad_or = text_reader_->ReadLine();
  if (!pad_or.ok()) {
    goto data_loss;
  }

  quality_or = text_reader_->ReadLine();
  if (!quality_or.ok()) {
    goto data_loss;
  }

  *header = header_or.ValueOrDie();
  *sequence = sequence_or.ValueOrDie();
  *pad = pad_or.ValueOrDie();
  *quality = quality_or.ValueOrDie();
  return ::nucleus::Status();

data_loss:
  return ::nucleus::DataLoss("Failed to parse FASTQ record");
}

StatusOr<std::shared_ptr<FastqIterable>> FastqReader::Iterate() const {
  if (!text_reader_) {
    return ::nucleus::FailedPrecondition(
        "Cannot Iterate a closed FastqReader.");
  }
  return StatusOr<std::shared_ptr<FastqIterable>>(
      MakeIterable<FastqFullFileIterable>(this));
}

// Iterable class definitions.
StatusOr<bool> FastqFullFileIterable::Next(FastqRecord* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  const FastqReader* fastq_reader = static_cast<const FastqReader*>(reader_);
  string header, sequence, pad, quality;
  ::nucleus::Status status =
      fastq_reader->Next(&header, &sequence, &pad, &quality);
  if (!status.ok()) {
    if (::nucleus::IsOutOfRange(status)) {
      return false;
    } else {
      return status;
    }
  }
  NUCLEUS_RETURN_IF_ERROR(ConvertToPb(header, sequence, pad, quality, out));
  return true;
}

FastqFullFileIterable::~FastqFullFileIterable() {}

FastqFullFileIterable::FastqFullFileIterable(const FastqReader* reader)
    : Iterable(reader) {}

}  // namespace nucleus
