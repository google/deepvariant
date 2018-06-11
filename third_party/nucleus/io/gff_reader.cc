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

// Implementation of gff_reader.h
#include "third_party/nucleus/io/gff_reader.h"

#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/strip.h"

#include "third_party/nucleus/protos/range.pb.h"

#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/logging.h"

constexpr char GFF_COMMENT_PREFIX[] = "#";

namespace tf = tensorflow;

using nucleus::genomics::v1::GffHeader;
using nucleus::genomics::v1::GffReaderOptions;
using nucleus::genomics::v1::GffRecord;

namespace nucleus {

namespace {

tf::Status ParseGffHeaderLine(const string& line, GffHeader* header) {
  if (absl::StartsWith(line, "##gff-version")) {
    // redacted
    // conversions
    header->set_gff_version(string(absl::StripPrefix(line, "##")));
  } else if (absl::StartsWith(line, "##sequence-region")) {
    std::vector<string> tokens = absl::StrSplit(line, ' ');
    if (tokens.size() != 4) {
      return tf::errors::DataLoss("Invalid sequence-region GFF header.");
    }
    // Parse seqid.
    string seqid = tokens[1];
    // Parse start, end.
    int64 start1, end1;
    if (!absl::SimpleAtoi(tokens[2], &start1)) {
      return tf::errors::Unknown("Can't parse GFF sequence-region start");
    }
    if (!absl::SimpleAtoi(tokens[3], &end1)) {
      return tf::errors::Unknown("Can't parse GFF sequence-region end");
    }
    // Convert to zero-based end-exclusive coordinates.
    int64 start = start1 - 1;
    int64 end = end1;
    // Write on header record.
    auto* sequence_region = header->add_sequence_regions();
    sequence_region->set_reference_name(seqid);
    sequence_region->set_start(start);
    sequence_region->set_end(end);
  } else {
    // redacted
  }

  return tf::Status::OK();
}

// Peeks into the GFF file to extract the header.
// redacted
tf::Status ReadGffHeader(const string& path, GffHeader* header) {
  CHECK(header != nullptr);
  header->Clear();

  StatusOr<std::unique_ptr<TextReader>> reader_or = TextReader::FromFile(path);
  TF_RETURN_IF_ERROR(reader_or.status());
  std::unique_ptr<TextReader> text_reader = std::move(reader_or.ValueOrDie());

  StatusOr<string> line_or;
  string line;
  while ((line_or = text_reader->ReadLine()).ok() &&
         absl::StartsWith(line = line_or.ValueOrDie(), GFF_COMMENT_PREFIX)) {
    TF_RETURN_IF_ERROR(ParseGffHeaderLine(line, header));
  }

  // Propagate error, if any.
  tf::Status status = line_or.status();
  if (!status.ok() && !tf::errors::IsOutOfRange(status)) {
    return status;
  }

  return tf::Status::OK();
}

tf::Status NextNonCommentLine(TextReader& text_reader, string* line) {
  CHECK(line != nullptr);
  string tmp;
  do {
    StatusOr<string> line_or = text_reader.ReadLine();
    TF_RETURN_IF_ERROR(line_or.status());
    tmp = line_or.ValueOrDie();
  } while (absl::StartsWith(tmp, GFF_COMMENT_PREFIX));

  *line = tmp;
  return tf::Status::OK();
}

tf::Status ConvertToPb(const string& line, GffRecord* record) {
  CHECK(record != nullptr);

  std::vector<string> fields = absl::StrSplit(line, '\t');
  if (fields.size() != 9) {
    return tf::errors::Unknown("Incorrect number of columns in a GFF record.");
  }

  // Parse line.
  string seq_id = fields[0];
  string source = fields[1];
  string type = fields[2];

  int64 start1, end1;
  if (!absl::SimpleAtoi(fields[3], &start1)) {
    return tf::errors::Unknown("Cannot parse GFF record `start`");
  }
  if (!absl::SimpleAtoi(fields[4], &end1)) {
    return tf::errors::Unknown("Cannot parse GFF record `end`");
  }
  // Convert to zero-based end-exclusive coordinate system.
  // redacted
  int64 start = start1 - 1;
  int64 end = end1;

  // Parse score.
  float score = 0.0f;
  if (fields[5] != ".") {
    if (!absl::SimpleAtof(fields[5], &score)) {
      return tf::errors::Unknown("Cannot parse GFF record `score`");
    }
  }
  // Parse strand.
  GffRecord::Strand strand;
  const string& strand_field = fields[6];
  if (strand_field == ".") {
    strand = GffRecord::UNSPECIFIED_STRAND;
  } else if (strand_field == "+") {
    strand = GffRecord::FORWARD_STRAND;
  } else if (strand_field == "-") {
    strand = GffRecord::REVERSE_STRAND;
  } else {
    return tf::errors::Unknown("Invalid GFF record `strand` encoding");
  }
  // Parse phase.
  // redacted
  int phase;
  const string& phase_field = fields[7];
  if (phase_field == ".") {
    phase = 0;
  } else if (!absl::SimpleAtoi(phase_field, &phase) ||
             !((phase >= 0) && (phase < 3))) {
    return tf::errors::Unknown("Invalid GFF record `phase` encoding.");
  }
  // redacted

  // Write on the record.
  record->Clear();
  record->mutable_range()->set_reference_name(seq_id);
  record->mutable_range()->set_start(start);
  record->mutable_range()->set_end(end);
  record->set_source(source);
  record->set_type(type);
  record->set_score(score);
  record->set_strand(strand);
  record->set_phase(phase);

  return tf::Status::OK();
}

}  // namespace

// ------- GFF iteration class

// Iterable class for traversing all GFF records in the file.
class GffFullFileIterable : public GffIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::GffRecord* out) override;

  // Constructor is invoked via GffReader::Iterate.
  GffFullFileIterable(const GffReader* reader);
  ~GffFullFileIterable() override;
};

// Iterable class definitions.
StatusOr<bool> GffFullFileIterable::Next(
    nucleus::genomics::v1::GffRecord* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  const GffReader* gff_reader = static_cast<const GffReader*>(reader_);
  string line;
  tf::Status status = NextNonCommentLine(*gff_reader->text_reader_, &line);
  if (tf::errors::IsOutOfRange(status)) {
    return false;
  } else {
    TF_RETURN_IF_ERROR(status);
  }
  TF_RETURN_IF_ERROR(ConvertToPb(line, out));
  return true;
}

GffFullFileIterable::~GffFullFileIterable() {}

GffFullFileIterable::GffFullFileIterable(const GffReader* reader)
    : Iterable(reader) {}

// ------- GFF reader class

StatusOr<std::unique_ptr<GffReader>> GffReader::FromFile(
    const string& gff_path,
    const nucleus::genomics::v1::GffReaderOptions& options) {
  StatusOr<std::unique_ptr<TextReader>> text_reader_or =
      TextReader::FromFile(gff_path);
  TF_RETURN_IF_ERROR(text_reader_or.status());

  GffHeader header;
  TF_RETURN_IF_ERROR(ReadGffHeader(gff_path, &header));

  return std::unique_ptr<GffReader>(
      new GffReader(std::move(text_reader_or.ValueOrDie()), options, header));
}

GffReader::GffReader(std::unique_ptr<TextReader> text_reader,
                     const GffReaderOptions& options, const GffHeader& header)
    : text_reader_(std::move(text_reader)),
      options_(options),
      header_(header) {}

StatusOr<std::shared_ptr<GffIterable>> GffReader::Iterate() const {
  if (!text_reader_)
    return tf::errors::FailedPrecondition("Cannot Iterate a closed GffReader.");
  return StatusOr<std::shared_ptr<GffIterable>>(
      MakeIterable<GffFullFileIterable>(this));
}

tensorflow::Status GffReader::Close() {
  if (!text_reader_) {
    return tf::errors::FailedPrecondition("GffReader already closed");
  }
  tf::Status status = text_reader_->Close();
  text_reader_ = nullptr;
  return status;
}

}  // namespace nucleus
