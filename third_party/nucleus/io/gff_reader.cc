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

// Implementation of gff_reader.h
#include "third_party/nucleus/io/gff_reader.h"

#include <limits>
#include <utility>
#include <vector>

#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/strings/strip.h"
#include "absl/types/optional.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/core/status.h"
#include "google/protobuf/map.h"

using nucleus::genomics::v1::GffHeader;
using nucleus::genomics::v1::GffReaderOptions;
using nucleus::genomics::v1::GffRecord;

namespace nucleus {

// Constants
// TODO: share these with the writer
constexpr char kGffCommentPrefix[] = "#";
constexpr char kGffMissingField[] = ".";
constexpr double kGffMissingDouble = -std::numeric_limits<double>::infinity();
constexpr int32 kGffMissingInt32 = -1;

namespace {

::nucleus::Status ParseGffHeaderLine(const string& line, GffHeader* header) {
  if (absl::StartsWith(line, "##gff-version")) {
    // TODO: get rid of pessimizing string_view -> string
    // conversions
    header->set_gff_version(string(absl::StripPrefix(line, "##")));
  } else if (absl::StartsWith(line, "##sequence-region")) {
    std::vector<string> tokens = absl::StrSplit(line, ' ');
    if (tokens.size() != 4) {
      return ::nucleus::DataLoss("Invalid sequence-region GFF header.");
    }
    // Parse seqid.
    string seqid = tokens[1];
    // Parse start, end.
    int64 start1, end1;
    if (!absl::SimpleAtoi(tokens[2], &start1)) {
      return ::nucleus::Unknown("Can't parse GFF sequence-region start");
    }
    if (!absl::SimpleAtoi(tokens[3], &end1)) {
      return ::nucleus::Unknown("Can't parse GFF sequence-region end");
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
    // TODO: other directives currently ignored.
  }

  return ::nucleus::Status();
}

// Peeks into the GFF file to extract the header.
// TODO: refactor header reading
::nucleus::Status ReadGffHeader(const string& path, GffHeader* header) {
  CHECK(header != nullptr);
  header->Clear();

  StatusOr<std::unique_ptr<TextReader>> reader_or = TextReader::FromFile(path);
  NUCLEUS_RETURN_IF_ERROR(reader_or.status());
  std::unique_ptr<TextReader> text_reader = std::move(reader_or.ValueOrDie());

  StatusOr<string> line_or;
  string line;
  while ((line_or = text_reader->ReadLine()).ok() &&
         absl::StartsWith(line = line_or.ValueOrDie(), kGffCommentPrefix)) {
    NUCLEUS_RETURN_IF_ERROR(ParseGffHeaderLine(line, header));
  }

  // Propagate error, if any.
  ::nucleus::Status status = line_or.status();
  if (!status.ok() && !::nucleus::IsOutOfRange(status)) {
    return status;
  }

  return ::nucleus::Status();
}

::nucleus::Status NextNonCommentLine(TextReader& text_reader, string* line) {
  CHECK(line != nullptr);
  string tmp;
  do {
    StatusOr<string> line_or = text_reader.ReadLine();
    NUCLEUS_RETURN_IF_ERROR(line_or.status());
    tmp = line_or.ValueOrDie();
  } while (absl::StartsWith(tmp, kGffCommentPrefix));

  *line = tmp;
  return ::nucleus::Status();
}

// Parses the text `attributes_string`, which is a ';'-delimited list
// of string-to-string '=' assignments, into a proto string->string map.
::nucleus::Status ParseGffAttributes(
    const string& attributes_string,
    google::protobuf::Map<string, string>* attributes_map) {
  if (attributes_string == kGffMissingField || attributes_string.empty()) {
    attributes_map->clear();
    return ::nucleus::Status();
  }

  google::protobuf::Map<string, string> tmp;
  std::vector<absl::string_view> assignments =
      absl::StrSplit(attributes_string, ';');
  for (absl::string_view assignment : assignments) {
    std::vector<absl::string_view> tokens = absl::StrSplit(assignment, '=');
    if (tokens.size() != 2) {
      return ::nucleus::Unknown("Cannot parse GFF attributes string");
    }
    // TODO: get rid of pessimizing string_view -> string
    // conversions
    string lhs = string(tokens[0]);
    string rhs = string(tokens[1]);
    tmp[lhs] = rhs;
  }
  *attributes_map = tmp;
  return ::nucleus::Status();
}

// Converts a text GFF line into a GffRecord proto message, or returns an error
// code if the line is malformed.  The record will only be modified if the call
// succeeds.
::nucleus::Status ConvertToPb(const string& line, GffRecord* record) {
  CHECK(record != nullptr);

  std::vector<string> fields = absl::StrSplit(line, '\t');
  if (fields.size() != 9) {
    return ::nucleus::Unknown("Incorrect number of columns in a GFF record.");
  }

  // Parse line.
  string seq_id = fields[0];
  if (seq_id == kGffMissingField || seq_id.empty()) {
    return ::nucleus::Unknown("GFF mandatory seq_id field is missing");
  }
  string source = (fields[1] == kGffMissingField ? "" : fields[1]);
  string type = (fields[2] == kGffMissingField ? "" : fields[2]);

  int64 start1, end1;
  if (!absl::SimpleAtoi(fields[3], &start1)) {
    return ::nucleus::Unknown("Cannot parse GFF record `start`");
  }
  if (!absl::SimpleAtoi(fields[4], &end1)) {
    return ::nucleus::Unknown("Cannot parse GFF record `end`");
  }
  // Convert to zero-based end-exclusive coordinate system.
  // TODO: consider factoring this out.
  int64 start = start1 - 1;
  int64 end = end1;

  // Parse score.
  absl::optional<float> score;
  if (fields[5] != kGffMissingField) {
    float value;
    if (!absl::SimpleAtof(fields[5], &value)) {
      return ::nucleus::Unknown("Cannot parse GFF record `score`");
    }
    score = value;
  }
  // Parse strand.
  GffRecord::Strand strand;
  const string& strand_field = fields[6];
  if (strand_field == kGffMissingField) {
    strand = GffRecord::UNSPECIFIED_STRAND;
  } else if (strand_field == "+") {
    strand = GffRecord::FORWARD_STRAND;
  } else if (strand_field == "-") {
    strand = GffRecord::REVERSE_STRAND;
  } else {
    return ::nucleus::Unknown("Invalid GFF record `strand` encoding");
  }
  // Parse phase.
  absl::optional<int> phase;
  const string& phase_field = fields[7];
  if (phase_field != kGffMissingField) {
    int value;
    if (!absl::SimpleAtoi(phase_field, &value) ||
        !((value >= 0) && (value < 3))) {
      return ::nucleus::Unknown("Invalid GFF record `phase` encoding.");
    }
    phase = value;
  }
  // Parse attributes dictionary.
  google::protobuf::Map<string, string> attributes_map;
  NUCLEUS_RETURN_IF_ERROR(ParseGffAttributes(fields[8], &attributes_map));

  // Write on the record.
  record->Clear();
  record->mutable_range()->set_reference_name(seq_id);
  record->mutable_range()->set_start(start);
  record->mutable_range()->set_end(end);
  record->set_source(source);
  record->set_type(type);
  record->set_score(score.value_or(kGffMissingDouble));
  record->set_strand(strand);
  record->set_phase(phase.value_or(kGffMissingInt32));
  *record->mutable_attributes() = attributes_map;

  return ::nucleus::Status();
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
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  const GffReader* gff_reader = static_cast<const GffReader*>(reader_);
  string line;
  ::nucleus::Status status =
      NextNonCommentLine(*gff_reader->text_reader_, &line);
  if (::nucleus::IsOutOfRange(status)) {
    return false;
  } else {
    NUCLEUS_RETURN_IF_ERROR(status);
  }
  NUCLEUS_RETURN_IF_ERROR(ConvertToPb(line, out));
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
  NUCLEUS_RETURN_IF_ERROR(text_reader_or.status());

  GffHeader header;
  NUCLEUS_RETURN_IF_ERROR(ReadGffHeader(gff_path, &header));

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
    return ::nucleus::FailedPrecondition("Cannot Iterate a closed GffReader.");
  return StatusOr<std::shared_ptr<GffIterable>>(
      MakeIterable<GffFullFileIterable>(this));
}

::nucleus::Status GffReader::Close() {
  if (!text_reader_) {
    return ::nucleus::FailedPrecondition("GffReader already closed");
  }
  ::nucleus::Status status = text_reader_->Close();
  text_reader_ = nullptr;
  return status;
}

}  // namespace nucleus
