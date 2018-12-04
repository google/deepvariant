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

#include "third_party/nucleus/io/unindexed_fasta_reader.h"

#include <stddef.h>
#include <utility>

#include "absl/strings/ascii.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/logging.h"

namespace tf = tensorflow;

namespace nucleus {

using genomics::v1::Range;

namespace {

// Helper method to get the name in a header line. This function assumes the
// first character is '>'.
absl::string_view GetNameInHeaderLine(absl::string_view line) {
  DCHECK_LT(1, line.size()) << "name must contain more than >";
  size_t space_idx = line.find(' ');
  if (space_idx == string::npos) {
    // No space is found. The name is the entire string after >.
    space_idx = line.size();
  }
  return line.substr(1, space_idx - 1);
}

}  // namespace

// Iterable class for traversing all Fasta records in the file.
class UnindexedFastaReaderIterable : public GenomeReferenceRecordIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(GenomeReferenceRecord* out) override;

  // Constructor is invoked via UnindexedFastaReader::Iterate.
  UnindexedFastaReaderIterable(const UnindexedFastaReader* reader);
  ~UnindexedFastaReaderIterable() override;

 private:
  // If non-empty, contains the name/id in the header line of the next record.
  std::string next_name_;
};

StatusOr<std::unique_ptr<UnindexedFastaReader>> UnindexedFastaReader::FromFile(
    const string& fasta_path) {
  StatusOr<std::unique_ptr<TextReader>> textreader_or =
      TextReader::FromFile(fasta_path);
  TF_RETURN_IF_ERROR(textreader_or.status());
  return std::unique_ptr<UnindexedFastaReader>(
      new UnindexedFastaReader(std::move(textreader_or.ValueOrDie())));
}

UnindexedFastaReader::~UnindexedFastaReader() {}

const std::vector<nucleus::genomics::v1::ContigInfo>&
UnindexedFastaReader::Contigs() const {
  LOG(FATAL) << "Unimplemented function invoked : " << __func__;
  return contigs_;
}

StatusOr<string> UnindexedFastaReader::GetBases(const Range& range) const {
  LOG(FATAL) << "Unimplemented function invoked : " << __func__;
  return tf::errors::Unimplemented("");
}

StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>
UnindexedFastaReader::Iterate() const {
  return StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>(
      MakeIterable<UnindexedFastaReaderIterable>(this));
}

tf::Status UnindexedFastaReader::Close() {
  if (!text_reader_) {
    return tf::errors::FailedPrecondition(
        "UnindexedFastaReader already closed");
  }
  // Close the file pointer.
  tf::Status close_status = text_reader_->Close();
  text_reader_ = nullptr;
  return close_status;
}

UnindexedFastaReader::UnindexedFastaReader(
    std::unique_ptr<TextReader> text_reader)
    : text_reader_(std::move(text_reader)) {}

StatusOr<bool> UnindexedFastaReaderIterable::Next(GenomeReferenceRecord* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  DCHECK(out && out->first.empty() && out->second.empty())
      << "out must be default initialized";

  const UnindexedFastaReader* fasta_reader =
      static_cast<const UnindexedFastaReader*>(reader_);
  if (!fasta_reader->text_reader_) {
    return tf::errors::FailedPrecondition(
        "Cannot iterate a closed UnindexedFastaReader.");
  }
  if (!next_name_.empty()) {
    out->first = next_name_;
    next_name_.clear();
  }
  bool eof = false;
  while (true) {
    // Read one line.
    StatusOr<string> line = fasta_reader->text_reader_->ReadLine();
    if (!line.ok()) {
      if (tf::errors::IsOutOfRange(line.status())) {
        eof = true;
        break;
      }
      return tf::errors::DataLoss("Failed to parse FASTA");
    }
    std::string l = line.ValueOrDie();

    if (l.empty()) continue;
    // Check if the line is a header or a sequence.
    if (l.at(0) == '>') {
      absl::string_view parsed_name = GetNameInHeaderLine(l);
      if (out->first.empty()) {
        out->first = string(parsed_name);
        continue;
      }
      next_name_ = string(parsed_name);
      return true;
    }
    // Processing a sequence line. If name is absent by now, return an error.
    if (out->first.empty()) {
      return tf::errors::DataLoss("Name not found in FASTA");
    }
    out->second.append(
        absl::AsciiStrToUpper(absl::StripTrailingAsciiWhitespace(l)));
  }
  if (eof && out->first.empty()) {
    // No more records.
    return false;
  }
  return true;
}

UnindexedFastaReaderIterable::~UnindexedFastaReaderIterable() {}

UnindexedFastaReaderIterable::UnindexedFastaReaderIterable(
    const UnindexedFastaReader* reader)
    : Iterable(reader) {}

}  // namespace nucleus
