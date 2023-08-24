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

#include "third_party/nucleus/io/sam_writer.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <memory>
#include <utility>

#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/strings/string_view.h"
#include "htslib/cram.h"
#include "htslib/hts_endian.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/sam_utils.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/core/status.h"
#include "google/protobuf/repeated_field.h"

namespace nucleus {

using genomics::v1::Read;
using genomics::v1::SamHeader;
using genomics::v1::Value;

namespace {

// Helper class to calculate byte representation (specified in SAM format) of
// the auxiliary info field of a Read message.
// Note that only string, int, and double fields are currently supported.
class AuxBuilder {
 public:
  // Type tags:
  // Null terminated string.
  const char kZTag = 'Z';
  // Signed int32.
  const char kiTag = 'i';
  // Null terminated string.
  const char kfTag = 'f';

  AuxBuilder(const Read& read) : read_(read) {}
  // Disable assignment/copy operations
  AuxBuilder(const AuxBuilder& other) = delete;
  AuxBuilder& operator=(const AuxBuilder&) = delete;

  // Returns the number of bytes needed to represent info fields in SAM format.
  StatusOr<size_t> NumBytes() {
    if (has_num_bytes_) {
      return num_bytes_;
    }
    auto status = ComputeNumBytes(read_);
    if (status.ok()) {
      has_num_bytes_ = true;
      num_bytes_ = status.ValueOrDie();
    }
    return status;
  }

  // Copies out the byte representation to |data|. This can only be called after
  // NumBytes() is successful and |data| must point to at least NumBytes()
  // bytes.
  void CopyTo(uint8_t* data) {
    // It's a programming error if CopyTo() is called not after a successful
    // NumBytes().
    CHECK(has_num_bytes_);
    uint8_t* data_array_ptr = data;
    for (const auto& entry : read_.info()) {
      if (entry.second.values_size() != 1) {
        // TODO: Support writing byte-array field.
        LOG(WARNING) << "SamWriter currently doesn't support writing info "
                        "fields of size "
                     << entry.second.values_size();
        continue;
      }
      const Value& v = entry.second.values(0);
      if (!(v.kind_case() == Value::kIntValue ||
            v.kind_case() == Value::kStringValue ||
            v.kind_case() == Value::kNumberValue)) {
        LOG(WARNING) << "unrecognized kind case " << v.kind_case();
        continue;
      }
      memcpy(data_array_ptr, entry.first.data(), 2);
      data_array_ptr += 2;
      switch (v.kind_case()) {
        case Value::kIntValue: {
          // Write as a 'i' (signed 32bit integer).
          memcpy(data_array_ptr, &kiTag, 1);
          i32_to_le(v.int_value(), data_array_ptr + 1);  // 4 bytes
          data_array_ptr += 5;
          break;
        }
        case Value::kStringValue: {
          // Write as a 'Z' (null-terminated string).
          memcpy(data_array_ptr, &kZTag, 1);
          const size_t string_size = v.string_value().size();
          memcpy(data_array_ptr + 1, v.string_value().data(), string_size);
          memset(data_array_ptr + 1 + string_size, 0, 1);
          // Leave one for null byte.
          data_array_ptr += 1 + string_size + 1;
          break;
        }
        case Value::kNumberValue: {
          // Write as a 'f' (4-byte float).
          memcpy(data_array_ptr, &kfTag, 1);
          double_to_le(v.number_value(), data_array_ptr + 1);
          data_array_ptr += 5;
          break;
        }
        default: {
          break;
        }
      }
    }
    DCHECK_EQ(data_array_ptr - data, num_bytes_);
  }

 private:
  static StatusOr<size_t> ComputeNumBytes(const Read& read) {
    size_t num_bytes = 0;
    for (const auto& entry : read.info()) {
      if (entry.first.size() != 2) {
        return ::nucleus::Unknown(absl::StrCat(
            "info key should be of two characters: ", entry.first));
      }
      if (entry.second.values_size() != 1) {
        // TODO: Support writing byte-array field.
        LOG(WARNING) << "SamWriter currently doesn't support writing info "
                        "fields of size "
                     << entry.second.values_size();
        continue;
      }
      const Value& v = entry.second.values(0);
      if (!(v.kind_case() == Value::kIntValue ||
            v.kind_case() == Value::kStringValue ||
            v.kind_case() == Value::kNumberValue)) {
        LOG(WARNING) << "unrecognized kind case " << v.kind_case();
        continue;
      }
      num_bytes += 3;
      // Two-character TAG and a one-character TYPE.
      switch (v.kind_case()) {
        case Value::kIntValue:
          // All integer types that are <= 32 bits (cCsSIi) were converted to
          // signed int32. internal.
          num_bytes += 4;
          break;
        case Value::kStringValue:
          num_bytes += entry.second.values(0).string_value().size() + 1;
          break;
        case Value::kNumberValue:
          // double to 4-byte float
          num_bytes += 4;
          break;
        default:
          break;
      }
    }
    return num_bytes;
  }
  const Read& read_;
  size_t num_bytes_ = 0;
  // Whether |num_bytes_| has been successfully computed and |read_| has been
  // validated through NumBytes().
  bool has_num_bytes_ = false;
};

// Helper method to append "\t{tag}{value}" to |text_stream|. |tag| should end
// with a colon.
void AppendTag(const char* tag, absl::string_view value, string* text_stream) {
  if (value.empty()) {
    return;
  }
  absl::StrAppend(text_stream, "\t", tag, value);
}

// Adds @HD lines to |text_stream| based on information in |sam_header|.
void AppendHeaderLine(const SamHeader& sam_header, string* text_stream) {
  absl::StrAppend(text_stream, kSamHeaderTag);
  AppendTag(kVNTag, sam_header.format_version(), text_stream);
  switch (sam_header.sorting_order()) {
    case SamHeader::UNKNOWN:
      AppendTag(kSOTag, "unknown", text_stream);
      break;
    case SamHeader::UNSORTED:
      AppendTag(kSOTag, "unsorted", text_stream);
      break;
    case SamHeader::QUERYNAME:
      AppendTag(kSOTag, "queryname", text_stream);
      break;
    case SamHeader::COORDINATE:
      AppendTag(kSOTag, "coordinate", text_stream);
      break;
    default:
      LOG(WARNING) << "unrecognized sorting order";
      break;
  }
  switch (sam_header.alignment_grouping()) {
    case SamHeader::NONE:
      AppendTag(kGOTag, "none", text_stream);
      break;
    case SamHeader::QUERY:
      AppendTag(kGOTag, "query", text_stream);
      break;
    case SamHeader::REFERENCE:
      AppendTag(kGOTag, "reference", text_stream);
      break;
    default:
      LOG(WARNING) << "unrecognized alignment group";
      break;
  }
  absl::StrAppend(text_stream, "\n");
}

// Adds @SQ lines to |text_stream| based on information in |sam_header|.
// Note that this is only required to be present for CRAM.  htslib's SAM/BAM
// header will write out @SQ lines automatically based on bam_hdr_t's
// target_name and target_len.
void AppendReferenceSequence(const SamHeader& sam_header, string* text_stream) {
  for (const auto& contig : sam_header.contigs()) {
    absl::StrAppend(text_stream, kSamReferenceSequenceTag);
    AppendTag(kSNTag, contig.name(), text_stream);
    AppendTag(kLNTag, std::to_string(contig.n_bases()), text_stream);
    absl::StrAppend(text_stream, "\n");
  }
}

// Adds @RG lines to |text_stream| based on information in |sam_header|.
void AppendReadGroups(const genomics::v1::SamHeader& sam_header,
                      string* text_stream) {
  for (const auto& read_group : sam_header.read_groups()) {
    absl::StrAppend(text_stream, kSamReadGroupTag);
    AppendTag(kIDTag, read_group.name(), text_stream);
    AppendTag(kCNTag, read_group.sequencing_center(), text_stream);
    AppendTag(kDSTag, read_group.description(), text_stream);
    AppendTag(kDTTag, read_group.date(), text_stream);
    AppendTag(kFOTag, read_group.flow_order(), text_stream);
    AppendTag(kKSTag, read_group.key_sequence(), text_stream);
    AppendTag(kLBTag, read_group.library_id(), text_stream);
    AppendTag(kPGTag,
              absl::StrJoin(read_group.program_ids().cbegin(),
                            read_group.program_ids().cend(), ","),
              text_stream);
    if (read_group.predicted_insert_size() != 0) {
      AppendTag(kPITag, std::to_string(read_group.predicted_insert_size()),
                text_stream);
    }
    AppendTag(kPLTag, read_group.platform(), text_stream);
    AppendTag(kPMTag, read_group.platform_model(), text_stream);
    AppendTag(kPUTag, read_group.platform_unit(), text_stream);
    AppendTag(kSMTag, read_group.sample_id(), text_stream);
    absl::StrAppend(text_stream, "\n");
  }
}

// Adds @PG lines to |text_stream| based on information in |sam_header|.
void AppendPrograms(const genomics::v1::SamHeader& sam_header,
                    string* text_stream) {
  for (const auto& program : sam_header.programs()) {
    absl::StrAppend(text_stream, kSamProgramTag);
    AppendTag(kIDTag, program.id(), text_stream);
    AppendTag(kPNTag, program.name(), text_stream);
    AppendTag(kCLTag, program.command_line(), text_stream);
    AppendTag(kPPTag, program.prev_program_id(), text_stream);
    AppendTag(kDSTag, program.description(), text_stream);
    AppendTag(kVNTag, program.version(), text_stream);
    absl::StrAppend(text_stream, "\n");
  }
}

// Appends @CO lines.
void AppendComments(const genomics::v1::SamHeader& sam_header,
                    string* text_stream) {
  for (const auto& comment : sam_header.comments()) {
    if (comment.empty()) continue;
    absl::StrAppend(text_stream, kSamCommentTag, "\t", comment, "\n");
  }
}

// Populates the fields in |h| based on information in |sam_header| proto.
::nucleus::Status PopulateNativeHeader(
    const genomics::v1::SamHeader& sam_header, bool is_cram, bam_hdr_t* h) {
  DCHECK_NE(nullptr, h);
  const uint32_t contig_size = sam_header.contigs().size();
  h->n_targets = contig_size;
  h->target_name = (char**)malloc(contig_size * sizeof(char*));
  h->target_len = (uint32_t*)malloc(contig_size * sizeof(uint32_t));

  int contig_i = 0;
  for (const auto& contig : sam_header.contigs()) {
    h->target_name[contig_i] =
        (char*)malloc((contig.name().size() + 1) * sizeof(char));
    memcpy(h->target_name[contig_i], contig.name().c_str(),
           contig.name().size());
    h->target_name[contig_i][contig.name().size()] = '\0';
    h->target_len[contig_i] = contig.n_bases();
    contig_i++;
  }

  // Consider avoiding the intermediate copy.
  string text;
  AppendHeaderLine(sam_header, &text);
  if (is_cram) {
    AppendReferenceSequence(sam_header, &text);
  }
  AppendReadGroups(sam_header, &text);
  AppendPrograms(sam_header, &text);
  AppendComments(sam_header, &text);

  h->l_text = text.length();
  h->text = (char*)malloc((text.length() + 1) * sizeof(char));
  memcpy(h->text, text.c_str(), text.length());
  h->text[text.length()] = '\0';

  return ::nucleus::Status();
}

// Returns a bam1_core_t.flag based on the information contained in |read|.
uint16_t GetReadFlag(const Read& read) {
  uint16_t flag = 0;
  if (read.proper_placement()) {
    flag |= BAM_FPROPER_PAIR;
  }
  if (read.duplicate_fragment()) {
    flag |= BAM_FDUP;
  }
  if (read.failed_vendor_quality_checks()) {
    flag |= BAM_FQCFAIL;
  }
  if (read.secondary_alignment()) {
    flag |= BAM_FSECONDARY;
  }
  if (read.supplementary_alignment()) {
    flag |= BAM_FSUPPLEMENTARY;
  }
  if (read.number_reads() == 2) {
    flag |= BAM_FPAIRED;
    if (read.read_number() == 0) {
      flag |= BAM_FREAD1;
    }
  }
  if (!read.has_next_mate_position()) {
    flag |= BAM_FMUNMAP;
  }
  if (!read.has_alignment()) {
    flag |= BAM_FUNMAP;
  }
  if (read.has_alignment() && read.alignment().position().reverse_strand()) {
    flag |= BAM_FREVERSE;
  }
  if (read.has_next_mate_position() &&
      read.next_mate_position().reverse_strand()) {
    flag |= BAM_FMREVERSE;
  }
  return flag;
}

// Populates the fields in |b| based on information in |h| and |read| proto.
::nucleus::Status PopulateNativeBody(const Read& read, const bam_hdr_t* h,
                                     bam1_t* b) {
  DCHECK_NE(nullptr, b);
  bam1_core_t* c = &b->core;
  c->isize = read.fragment_length();
  c->flag = GetReadFlag(read);
  c->l_qseq = read.aligned_sequence().size();
  if (read.has_alignment()) {
    c->qual = read.alignment().mapping_quality();
    if (read.alignment().has_position()) {
      c->pos = read.alignment().position().position();
      for (int i = 0; i < h->n_targets; ++i) {
        if (h->target_name[i] == read.alignment().position().reference_name()) {
          c->tid = i;
          break;
        }
      }
    }
  }
  if (read.has_next_mate_position()) {
    c->mpos = read.next_mate_position().position();
    if (read.next_mate_position().reference_name() == "*") {
      c->mtid = -1;
    } else {
      for (int i = 1; i < h->n_targets; ++i) {
        if (h->target_name[i] == read.next_mate_position().reference_name()) {
          c->mtid = i;
          break;
        }
      }
    }
  }
  // |b->l_data| is the length of concatenated structure:
  // qname-cigar-seq-qual-aux.

  // Copy qname with a null terminator.
  c->l_qname = read.fragment_name().size() + 1;
  // Each cigar takes 4 bytes.
  const size_t cigar_bytes =
      4 * (read.has_alignment() ? read.alignment().cigar_size() : 0);
  // Each base is represented using 4 bit, so one byte can represent 2 bases.
  const size_t encoded_base_bytes = (read.aligned_sequence().size() + 1) >> 1;
  // Each qual is 1 byte.
  const size_t aligned_quality_bytes = read.aligned_quality_size();
  // Use a helper class to calculate the number of bytes of all auxiliary info.
  AuxBuilder auxBuilder(read);

  StatusOr<size_t> aux_status = auxBuilder.NumBytes();
  size_t aux_bytes = 0;
  if (aux_status.ok()) {
    aux_bytes = aux_status.ValueOrDie();
  }
  size_t data_array_bytes = c->l_qname + cigar_bytes + encoded_base_bytes +
                            aligned_quality_bytes + aux_bytes;
  if (read.has_alignment()) {
    c->n_cigar = read.alignment().cigar_size();
  }

  // array is freed by htslib
  auto data_array = (uint8_t*)malloc(data_array_bytes);
  b->l_data = data_array_bytes;
  b->m_data = data_array_bytes;
  b->data = data_array;

  auto data_array_ptr = data_array;

  // Copy qname.
  memcpy(data_array_ptr, read.fragment_name().c_str(), c->l_qname);
  data_array_ptr += c->l_qname;

  // Copy cigars.
  for (const auto& cigar : read.alignment().cigar()) {
    uint32_t value =
        (cigar.operation_length() << BAM_CIGAR_SHIFT) |
        kProtoToHtslibCigar[static_cast<uint32_t>(cigar.operation())];
    memcpy(data_array_ptr, &value, 4);
    data_array_ptr += 4;
  }

  // Copy bases.
  memset(data_array_ptr, 0, encoded_base_bytes);
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    // seq_nt16_table is defined in htslib/hts.h and gives us the 8 bit encoded
    // value of each base.
    uint8_t encoded = seq_nt16_table[(int)read.aligned_sequence()[i]];
    data_array_ptr[i >> 1] |= encoded << ((~i & 1) << 2);
  }
  data_array_ptr += encoded_base_bytes;

  // Copy qual.
  for (const auto& qual : read.aligned_quality()) {
    memcpy(data_array_ptr, &qual, 1);
    data_array_ptr += 1;
  }

  if (aux_status.ok()) {
    auxBuilder.CopyTo(data_array_ptr);
  }
  return ::nucleus::Status();
}

// Helper method to get file extension of |file_path|.
string GetFileExtension(absl::string_view file_path) {
  auto pos = file_path.rfind('.');
  if (pos == absl::string_view::npos) {
    return string();
  }
  return string(file_path.substr(pos + 1));
}

}  // namespace

// -----------------------------------------------------------------------------
//
// Wrapper classes to take ownership of native htslib objects.
//
// -----------------------------------------------------------------------------

class SamWriter::NativeHeader {
 public:
  NativeHeader(bam_hdr_t* h) : h_(h) {}
  ~NativeHeader() { bam_hdr_destroy(h_); }
  // Disable assignment/copy operations
  NativeHeader(const NativeHeader& other) = delete;
  NativeHeader& operator=(const NativeHeader&) = delete;

  bam_hdr_t* value() { return h_; }

 private:
  bam_hdr_t* const h_;
};

class SamWriter::NativeFile {
 public:
  NativeFile(samFile* f) : f_(f) {}
  ~NativeFile() { hts_close(f_); }
  // Disable assignment/copy operations
  NativeFile(const NativeFile& other) = delete;
  NativeFile& operator=(const NativeFile&) = delete;

  samFile* value() { return f_; }

 private:
  samFile* const f_;
};

class SamWriter::NativeBody {
 public:
  NativeBody(bam1_t* b) : b_(b) {}
  ~NativeBody() { bam_destroy1(b_); }
  // Disable assignment/copy operations
  NativeBody(const NativeBody& other) = delete;
  NativeBody& operator=(const NativeBody&) = delete;

  bam1_t* value() { return b_; }

 private:
  bam1_t* const b_;
};

// -----------------------------------------------------------------------------
//
// Writer for SAM formats containing NGS reads.
//
// -----------------------------------------------------------------------------

StatusOr<std::unique_ptr<SamWriter>> SamWriter::ToFile(
    const string& sam_path, const genomics::v1::SamHeader& sam_header) {
  return ToFile(sam_path, string(), false, sam_header);
}

StatusOr<std::unique_ptr<SamWriter>> SamWriter::ToFile(
    const string& sam_path, const string& ref_path, bool embed_ref,
    const genomics::v1::SamHeader& sam_header) {
  htsFormat fmt;
  fmt.specific = nullptr;

  if (hts_parse_format(&fmt, GetFileExtension(sam_path).c_str()) < 0) {
    return ::nucleus::Unknown(
        absl::StrCat("Parsing file format fails: ", sam_path));
  }
  samFile* fp = hts_open_format_x(sam_path, "w", &fmt);
  if (fp == nullptr) {
    return ::nucleus::Unknown(
        absl::StrCat("Could not open file for writing: ", sam_path));
  }
  // Set user provided reference FASTA to decode CRAM.
  if (fp->format.format == cram) {
    if (ref_path.empty()) {
      return ::nucleus::FailedPrecondition(
          "Writing CRAM format requires a reference file");
    }
    LOG(INFO) << "Setting CRAM reference path to '" << ref_path << "'";
    if (cram_set_option(fp->fp.cram, CRAM_OPT_REFERENCE, ref_path.c_str()) <
        0) {
      return ::nucleus::Unknown(absl::StrCat(
          "Failed to set the CRAM_OPT_REFERENCE value to ", ref_path));
    }
    cram_set_option(fp->fp.cram, CRAM_OPT_EMBED_REF, embed_ref ? 1 : 0);
  }

  auto native_file = std::make_unique<NativeFile>(fp);
  auto native_header = std::make_unique<NativeHeader>(bam_hdr_init());
  NUCLEUS_RETURN_IF_ERROR(PopulateNativeHeader(
      sam_header, fp->format.format == cram, native_header->value()));

  if (sam_hdr_write(fp, native_header->value()) < 0) {
    return ::nucleus::Unknown("Writing header to file failed");
  }
  return absl::WrapUnique<SamWriter>(
      new SamWriter(std::move(native_file), std::move(native_header)));
}

SamWriter::SamWriter(std::unique_ptr<NativeFile> file,
                     std::unique_ptr<NativeHeader> header)
    : native_file_(std::move(file)), native_header_(std::move(header)) {}

SamWriter::~SamWriter() {
  if (native_file_) {
    // There's nothing we can do but assert fail if there's an error during
    // the Close() call here.
    NUCLEUS_CHECK_OK(Close());
  }
}

::nucleus::Status SamWriter::Close() {
  native_file_.reset();
  native_header_ = nullptr;
  return ::nucleus::Status();
}

::nucleus::Status SamWriter::Write(const Read& read) {
  auto body = std::make_unique<NativeBody>(bam_init1());
  ::nucleus::Status status =
      PopulateNativeBody(read, native_header_->value(), body->value());
  if (!status.ok()) {
    return status;
  }
  if (sam_write1(native_file_->value(), native_header_->value(),
                 body->value()) < 0) {
    return ::nucleus::Unknown("Cannot add record");
  }
  return ::nucleus::Status();
}

}  // namespace nucleus
