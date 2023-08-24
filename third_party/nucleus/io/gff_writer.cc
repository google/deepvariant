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

#include "third_party/nucleus/io/gff_writer.h"

#include <limits>
#include <map>
#include <utility>

#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/strings/substitute.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"
#include "google/protobuf/map.h"

namespace nucleus {

using genomics::v1::GffHeader;
using genomics::v1::GffRecord;
using genomics::v1::GffWriterOptions;
using genomics::v1::Range;

// Constants
// TODO: share these with the reader.
constexpr char kGffCommentPrefix[] = "#";
constexpr char kGffMissingField[] = ".";
constexpr double kGffMissingDouble = -std::numeric_limits<double>::infinity();
constexpr int32 kGffMissingInt32 = -1;

namespace {

::nucleus::Status WriteGffHeader(const GffHeader& header,
                                 TextWriter* text_writer) {
  NUCLEUS_RETURN_IF_ERROR(text_writer->Write("##gff-version 3.2.1\n"));
  for (const Range& range : header.sequence_regions()) {
    NUCLEUS_RETURN_IF_ERROR(text_writer->Write(
        // Range start converted from 0- to 1-based, end-inclusive.
        absl::Substitute("##sequence-region $0 $1 $2\n", range.reference_name(),
                         range.start() + 1, range.end())));
  }
  // TODO: write ontology headers.
  return ::nucleus::Status();
}

::nucleus::Status FormatGffAttributes(const GffRecord& record,
                                      string* gff_attributes) {
  // Sort to ensure deterministic iteration order.
  std::map<string, string> sorted_attributes(record.attributes().begin(),
                                             record.attributes().end());
  *gff_attributes =
      absl::StrJoin(sorted_attributes, ";", absl::PairFormatter("="));
  return ::nucleus::Status();
}

::nucleus::Status FormatGffLine(const GffRecord& record, string* gff_line) {
  string tmp, attributes;
  absl::StrAppend(&tmp, record.range().reference_name(), "\t");
  absl::StrAppend(
      &tmp, (record.source().empty() ? kGffMissingField : record.source()),
      "\t");
  absl::StrAppend(
      &tmp, (record.type().empty() ? kGffMissingField : record.type()), "\t");
  // Convert range to 1-based coordinates for GFF text.
  int64 start1 = record.range().start() + 1;
  int64 end1 = record.range().end();
  absl::StrAppend(&tmp, start1, "\t", end1, "\t");
  // Score
  string score_str =
      (record.score() == kGffMissingDouble ? kGffMissingField
                                           : absl::StrCat(record.score()));
  absl::StrAppend(&tmp, score_str, "\t");
  // Strand
  string strand_code;
  switch (record.strand()) {
    case GffRecord::UNSPECIFIED_STRAND:
      strand_code = kGffMissingField;
      break;
    case GffRecord::FORWARD_STRAND:
      strand_code = "+";
      break;
    case GffRecord::REVERSE_STRAND:
      strand_code = "-";
      break;
    default:
      return ::nucleus::InvalidArgument("Illegal GffRecord strand encoding");
  }
  absl::StrAppend(&tmp, strand_code, "\t");
  // Phase
  int phase = record.phase();
  if (phase >= 0 && phase <= 2) {
    absl::StrAppend(&tmp, phase, "\t");
  } else if (phase == kGffMissingInt32) {
    absl::StrAppend(&tmp, kGffMissingField, "\t");
  } else {
    return ::nucleus::InvalidArgument("Illegal GffRecord phase encoding");
  }
  // Attributes
  NUCLEUS_RETURN_IF_ERROR(FormatGffAttributes(record, &attributes));
  absl::StrAppend(&tmp, attributes);
  absl::StrAppend(&tmp, "\n");

  *gff_line = tmp;
  return ::nucleus::Status();
}

}  // namespace

StatusOr<std::unique_ptr<GffWriter>> GffWriter::ToFile(
    const string& gff_path, const GffHeader& header,
    const GffWriterOptions& options) {
  StatusOr<std::unique_ptr<TextWriter>> text_writer_or =
      TextWriter::ToFile(gff_path);
  NUCLEUS_RETURN_IF_ERROR(text_writer_or.status());

  std::unique_ptr<TextWriter> text_writer =
      std::move(text_writer_or.ValueOrDie());
  NUCLEUS_RETURN_IF_ERROR(WriteGffHeader(header, text_writer.get()));

  return absl::WrapUnique(
      new GffWriter(std::move(text_writer), header, options));
}

::nucleus::Status GffWriter::Write(const GffRecord& record) {
  if (!text_writer_)
    return ::nucleus::FailedPrecondition("Cannot write to closed GFF stream.");
  string line;
  NUCLEUS_RETURN_IF_ERROR(FormatGffLine(record, &line));
  return text_writer_->Write(line);
}

::nucleus::Status GffWriter::Close() {
  if (!text_writer_)
    return ::nucleus::FailedPrecondition(
        "Cannot close an already closed GffWriter");
  // Close the file pointer we have been writing to.
  ::nucleus::Status close_status = text_writer_->Close();
  text_writer_ = nullptr;
  return close_status;
}

GffWriter::GffWriter(std::unique_ptr<TextWriter> text_writer,
                     const GffHeader& header, const GffWriterOptions& options)
    : header_(header),
      options_(options),
      text_writer_(std::move(text_writer)) {}

}  // namespace nucleus
