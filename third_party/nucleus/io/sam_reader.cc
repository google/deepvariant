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

// Implementation of sam_reader.h
#include "third_party/nucleus/io/sam_reader.h"

#include <errno.h>
#include <stdint.h>

#include <algorithm>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "htslib/cram.h"
#include "htslib/hts.h"
#include "htslib/hts_endian.h"
#include "htslib/sam.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/sam_utils.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"
#include "google/protobuf/repeated_field.h"

namespace nucleus {

using absl::string_view;
using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::CigarUnit_Operation;
using nucleus::genomics::v1::Position;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::SamHeader;
using nucleus::genomics::v1::SamReaderOptions;
using std::vector;

using google::protobuf::RepeatedField;

namespace {

inline constexpr char kOQ[] = "OQ";

bool FileTypeIsIndexable(htsFormat format) {
  return format.format == bam || format.format == cram;
}

void AddHeaderLineToHeader(const string& line, SamHeader& header) {
  int tagLen = 3;

  static const std::map<string, nucleus::genomics::v1::SamHeader_SortingOrder>
      sorting_order_map = {{"coordinate", SamHeader::COORDINATE},
                           {"queryname", SamHeader::QUERYNAME},
                           {"unknown", SamHeader::UNKNOWN},
                           {"unsorted", SamHeader::UNSORTED}};

  static const std::map<string,
                        nucleus::genomics::v1::SamHeader_AlignmentGrouping>
      alignment_grouping_map = {{"none", SamHeader::NONE},
                                {"query", SamHeader::QUERY},
                                {"reference", SamHeader::REFERENCE}};

  for (const string_view& token : absl::StrSplit(line, '\t')) {
    if (token == kSamHeaderTag) continue;
    // TODO: remove string conversions.
    const string tag = string(token.substr(0, tagLen));
    const string value = string(token.substr(tagLen));
    if (tag == kVNTag) {
      header.set_format_version(value);
    } else if (tag == kSOTag) {
      const auto& it = sorting_order_map.find(value);
      if (it == sorting_order_map.end()) {
        LOG(WARNING) << "Unknown sorting order, defaulting to unknown: "
                     << line;
        header.set_sorting_order(SamHeader::UNKNOWN);
      } else {
        header.set_sorting_order(it->second);
      }
    } else if (tag == kGOTag) {
      const auto& it = alignment_grouping_map.find(value);
      if (it == alignment_grouping_map.end()) {
        LOG(WARNING) << "Unknown alignment grouping, defaulting to none: "
                     << line;
        header.set_alignment_grouping(SamHeader::NONE);
      } else {
        header.set_alignment_grouping(it->second);
      }
    } else if (tag == kSamHeaderTag) {
      // Skip this since it's the header line token.
    } else {
      LOG(WARNING) << "Unknown tag " << tag
                   << " in header line, ignoring: " << line;
    }
  }
}

void AddReadGroupToHeader(const string& line,
                          nucleus::genomics::v1::ReadGroup* readgroup) {
  int tagLen = 3;
  for (const string_view token : absl::StrSplit(line, '\t')) {
    if (token == kSamReadGroupTag) continue;
    // TODO: remove string conversions.
    const string tag = string(token.substr(0, tagLen));
    const string value = string(token.substr(tagLen));
    if (tag == kIDTag) {
      readgroup->set_name(value);
    } else if (tag == kCNTag) {
      readgroup->set_sequencing_center(value);
    } else if (tag == kDSTag) {
      readgroup->set_description(value);
    } else if (tag == kDTTag) {
      readgroup->set_date(value);
    } else if (tag == kFOTag) {
      readgroup->set_flow_order(value);
    } else if (tag == kKSTag) {
      readgroup->set_key_sequence(value);
    } else if (tag == kLBTag) {
      readgroup->set_library_id(value);
    } else if (tag == kPGTag) {
      readgroup->add_program_ids(value);
    } else if (tag == kPITag) {
      int size;
      CHECK(absl::SimpleAtoi(value, &size));
      readgroup->set_predicted_insert_size(size);
    } else if (tag == kPLTag) {
      readgroup->set_platform(value);
    } else if (tag == kPMTag) {
      readgroup->set_platform_model(value);
    } else if (tag == kPUTag) {
      readgroup->set_platform_unit(value);
    } else if (tag == kSMTag) {
      readgroup->set_sample_id(value);
    } else {
      LOG(WARNING) << "Unknown tag " << tag
                   << " in RG line, ignoring: " << line;
    }
  }
}

void AddProgramToHeader(const string& line,
                        nucleus::genomics::v1::Program* program) {
  int tagLen = 3;
  for (const string_view token : absl::StrSplit(line, '\t')) {
    if (token == kSamProgramTag) continue;
    // TODO: remove string conversions.
    const string tag = string(token.substr(0, tagLen));
    const string value = string(token.substr(tagLen));
    if (tag == kIDTag) {
      program->set_id(value);
    } else if (tag == kPNTag) {
      program->set_name(value);
    } else if (tag == kCLTag) {
      program->set_command_line(value);
    } else if (tag == kPPTag) {
      program->set_prev_program_id(value);
    } else if (tag == kDSTag) {
      program->set_description(value);
    } else if (tag == kVNTag) {
      program->set_version(value);
    }
  }
}

// This function is only checking the Read fields that has been filled.
// If you're modifying this function, please make sure any fields here have
// been previously filled.
bool PartialReadSatisfiesRequirements(
    const Read& read,
    const nucleus::genomics::v1::ReadRequirements& requirements) {
  // Using the fields that have been filled so far to check a subset of the
  // requirement in ReadSatisfiesRequirements function in utils.cc.
  // This helps us to abort earlier as needed.
  return (requirements.keep_duplicates() || !read.duplicate_fragment()) &&
         (requirements.keep_failed_vendor_quality_checks() ||
          !read.failed_vendor_quality_checks()) &&
         (requirements.keep_secondary_alignments() ||
          !read.secondary_alignment()) &&
         (requirements.keep_supplementary_alignments() ||
          !read.supplementary_alignment());
}

}  // namespace

namespace sam_reader_internal {
// Returns false if Read does not satisfy all of the ReadRequirements.
bool ReadSatisfiesRequirements(
    const Read& read,
    const nucleus::genomics::v1::ReadRequirements& requirements) {
  return PartialReadSatisfiesRequirements(read, requirements) &&
         (requirements.keep_unaligned() || read.has_alignment()) &&
         (requirements.keep_improperly_placed() ||
          IsReadProperlyPlaced(read)) &&
         (!read.has_alignment() || read.alignment().mapping_quality() >=
                                       requirements.min_mapping_quality());
}
}  // namespace sam_reader_internal

// -----------------------------------------------------------------------------
//
// Reader for SAM/BAM/CRAM etc formats containing NGS reads supported by htslib.
//
// -----------------------------------------------------------------------------

// Gets the size in bytes for a SAM/BAM aux tag based on their declared type.
//
// Based on code from htslib/sam.c which isn't exported from htslib. Returns
// a value < 0 if type isn't one of the expected types from htslib.
static inline int HtslibAuxSize(uint8_t type) {
  switch (type) {
    case 'A':
    case 'c':
    case 'C':
      return 1;
    case 's':
    case 'S':
      return 2;
    case 'f':
    case 'i':
    case 'I':
      return 4;
    default:
      return -1;  // -1 indicates error here.
  }
}

// Returns true iff query starts with the first prefix_len letters of prefix.
static inline bool StartsWith(const string& query, const char prefix[],
                              int prefix_len) {
  return query.compare(0, prefix_len, prefix) == 0;
}

// Parses out the aux tag attributes of a SAM record.
//
// From https://samtools.github.io/hts-specs/SAMv1.pdf
// 1.5 The alignment section: optional fields
//
// All optional fields follow the TAG:TYPE:VALUE format where TAG is a
// two-character string that matches /[A-Za-z][A-Za-z0-9]/. Each TAG can only
// appear once in one alignment line. A TAG containing lowercase letters are
// reserved for end users. In an optional field, TYPE is a single
// case-sensitive letter which defines the format of VALUE:
//
// Type Regexp matching VALUE Description
// A [!-~] Printable character
// i [-+]?[0-9]+ Signed integer
// f [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? Single-precision floating number
// Z [ !-~]* Printable string, including space
// H ([0-9A-F][0-9A-F])* Byte array in the Hex format
// B [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ Integer or numeric
// array
//
// For an integer or numeric array (type ‘B’), the first letter indicates the
// type of numbers in the following comma separated array. The letter can be
// one of ‘cCsSiIf’, corresponding to int8 t (signed 8-bit integer), uint8 t
// (unsigned 8-bit integer), int16 t, uint16 t, int32 t, uint32 t and float,
// respectively.
//
// Args:
//   b: The htslib bam record we will parse aux fields from.
//   option: Controls how aux fields are parsed.
//   read_message: Destination for parsed aux fields.
//
// Returns:
//   tensorflow::Status. Will be ok() if parsing succeeded or was not required,
//   otherwise will contain an error_message describing the problem.
::nucleus::Status ParseAuxFields(const bam1_t* b,
                                 const SamReaderOptions& options,
                                 Read* read_message) {
  if (options.aux_field_handling() != SamReaderOptions::PARSE_ALL_AUX_FIELDS) {
    return ::nucleus::Status();
  }
  const auto& aux_fields_to_keep = options.aux_fields_to_keep();

  uint8_t* s = bam_get_aux(b);
  const uint8_t* end = b->data + b->l_data;
  while (end - s >= 4) {
    // Each block is encoded like (each element is a byte):
    // [tag char 1, tag char 2, type byte, ...]
    // where the ... contents depends on the 2-character tag and type.
    const string tag = string(reinterpret_cast<char*>(s), 2);
    bool include_tag =
        (aux_fields_to_keep.empty() ||
         std::find(aux_fields_to_keep.begin(), aux_fields_to_keep.end(), tag) !=
             aux_fields_to_keep.end());
    s += 2;
    const uint8_t type = *s++;
    switch (type) {
      // An 'A' is just a single character string.
      case 'A': {
        if (include_tag) {
          // Safe since we know s is at least 4 bytes from the end.
          const string value = string(reinterpret_cast<char*>(s), 1);
          SetInfoField(tag, value, read_message);
        }
        s += 1;
      } break;
      // These are all different byte-sized integers.
      case 'C':
      case 'c':
      case 'S':
      case 's':
      case 'I':
      case 'i': {
        const int size = HtslibAuxSize(type);
        if (size < 0 || end - s < size)
          return ::nucleus::DataLoss("Malformed tag " + tag);
        errno = 0;
        const int value = bam_aux2i(s - 1);
        if (value == 0 && errno == EINVAL)
          return ::nucleus::DataLoss("Malformed tag " + tag);
        if (include_tag) SetInfoField(tag, value, read_message);
        s += size;
      } break;
      // A 4-byte floating point.
      case 'f': {
        if (end - s < 4) return ::nucleus::DataLoss("Malformed tag " + tag);
        const float value = le_to_float(s);
        if (include_tag) SetInfoField(tag, value, read_message);
        s += 4;
      } break;
      // Z and H are null-terminated strings.
      case 'Z':
      case 'H': {
        char* value = reinterpret_cast<char*>(s);
        for (; s < end && *s; ++s) {
        }  // Loop to the end.
        if (s >= end) return ::nucleus::DataLoss("Malformed tag " + tag);
        s++;
        // The H hex tag is not really used and likely deprecated (see:
        // https://sourceforge.net/p/samtools/mailman/message/28274509/
        // so we are explicitly skipping them here.
        if (type == 'Z' && include_tag) SetInfoField(tag, value, read_message);
      } break;
      // B is an array of atomic types (strings, ints, floats).
      case 'B': {
        const uint8_t sub_type = *s++;
        const int element_size = HtslibAuxSize(sub_type);
        if (element_size < 0)
          return ::nucleus::DataLoss("element_size == 0 for tag " + tag);
        // Prevents us from reading off the end of our buffer with le_to_u32.
        if (end - s < 4)
          return ::nucleus::DataLoss("data too short for tag " + tag);
        const int n_elements = le_to_u32(s);
        if (n_elements == 0) return ::nucleus::DataLoss("n_elements is zero");
        // We need to skip 4 bytes for n_elements int that occurs before the
        // array.
        s += 4;
        if (sub_type == 'c') {
          std::vector<int8_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              int8_t value = le_to_i8(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 'C') {
          std::vector<uint8_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              uint8_t value = reinterpret_cast<uint8_t>(*s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 's') {
          std::vector<int16_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              int16_t value = le_to_i16(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 'S') {
          std::vector<uint16_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              uint16_t value = le_to_u16(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 'i') {
          std::vector<int32_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              int32_t value = le_to_i32(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 'I') {
          std::vector<uint32_t> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              uint32_t value = le_to_u32(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else if (sub_type == 'f') {
          std::vector<float> all_values;
          for (int i = 0; i < n_elements; i++) {
            if (include_tag) {
              float value = le_to_float(s);
              all_values.push_back(value);
            }
            s += element_size;
          }
          if (include_tag) SetInfoField(tag, all_values, read_message);
        } else {
          return ::nucleus::DataLoss("Unknown subtype " +
                                     std::to_string(sub_type));
        }
      } break;
      default: {
        return ::nucleus::DataLoss("Unknown tag " + tag);
      }
    }
  }

  // Everything parsed correctly, so we return OK.
  return ::nucleus::Status();
}

// Assign aligned_quality. Depending on the use_original_base_quality_scores
// aligned_quality is read either from "QUAL" field or from "OQ" tag in SAM/BAM.
::nucleus::Status AssignAlignedQuality(const bam1_t* b,
                                       const SamReaderOptions& options,
                                       Read* read_message) {
  const bam1_core_t* c = &b->core;
  // Use optional "OQ" tag.
  if (options.use_original_base_quality_scores()) {
    const auto& info = read_message->info();
    auto info_it = info.find(kOQ);
    if (info_it != read_message->info().end() &&
        !info_it->second.values().empty()) {
      RepeatedField<int32>* quality = read_message->mutable_aligned_quality();
      quality->Reserve(c->l_qseq);
      const auto& oq_tag_value = *(info_it->second.values().begin());
      for (char c : oq_tag_value.string_value()) {
        quality->Add(reinterpret_cast<int>(c - 33));
      }
      return ::nucleus::Status();
    }
  } else {  // Use "QUAL" field.
    if (c->l_qseq) {
      uint8_t* quals = bam_get_qual(b);
      if (quals[0] != 0xff) {  // Not missing
        // TODO: Is there a more efficient way to do this?
        RepeatedField<int32>* quality = read_message->mutable_aligned_quality();
        quality->Reserve(c->l_qseq);
        for (int i = 0; i < c->l_qseq; ++i) {
          quality->Add(quals[i]);
        }
        return ::nucleus::Status();
      }
    }
  }
  return ::nucleus::Status(absl::StatusCode::kNotFound,
                           "Could not read base quality scores");
}

// Returns with tensorflow::error::Code::ABORTED if the read doesn't
// satisfy read requirements. When that happens, the function aborts early and
// doesn't fill the other fields such as aligned_sequence, which can be
// expensive in long reads.
::nucleus::Status ConvertToPb(const bam_hdr_t* h, const bam1_t* b,
                              const SamReaderOptions& options,
                              Read* read_message) {
  CHECK(h != nullptr) << "BAM header cannot be null";
  CHECK(b != nullptr) << "BAM record cannot be null";
  CHECK(read_message != nullptr) << "Read record cannot be null";

  read_message->Clear();

  const bam1_core_t* c = &b->core;

  // Grab a bunch of basic information from the bam1_t record and put it into
  // our protobuf.
  read_message->set_fragment_name(bam_get_qname(b));
  read_message->set_fragment_length(c->isize);
  read_message->set_proper_placement(c->flag & BAM_FPROPER_PAIR);
  read_message->set_duplicate_fragment(c->flag & BAM_FDUP);
  read_message->set_failed_vendor_quality_checks(c->flag & BAM_FQCFAIL);
  read_message->set_secondary_alignment(c->flag & BAM_FSECONDARY);
  read_message->set_supplementary_alignment(c->flag & BAM_FSUPPLEMENTARY);

  // Set the pairing information, which has a bit of complex logic to it. The
  // read number and number_reads per fragment depends on whether the read is
  // paired and, if so, whether its the first or second read.
  bool paired = c->flag & BAM_FPAIRED;
  read_message->set_read_number(c->flag & BAM_FREAD1 || !paired ? 0 : 1);
  read_message->set_number_reads(paired ? 2 : 1);

  if (options.has_read_requirements() &&
      !PartialReadSatisfiesRequirements(*read_message,
                                        options.read_requirements())) {
    return ::nucleus::Status(absl::StatusCode::kAborted,
                             "Read doesn't satisfy requirements.");
  }

  if (c->l_qseq) {
    // Convert the seq if it is present.
    string* read_seq = read_message->mutable_aligned_sequence();
    read_seq->reserve(c->l_qseq);
    uint8_t* seq = bam_get_seq(b);  // seq is stored as 8-bit offsets.
    for (int i = 0; i < c->l_qseq; ++i) {
      // Convert the offsets to upper case characters by their offset in
      // the constant seq_nt16_str from htslib.
      read_seq->push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
  }

  if (!(c->flag & BAM_FUNMAP)) {
    // If the read is mapped, set the mapping information in read_message.
    auto* linear_alignment = read_message->mutable_alignment();
    linear_alignment->set_mapping_quality(c->qual);

    if (c->n_cigar) {  // Convert our Cigar.
      uint32* cigar = bam_get_cigar(b);
      for (uint32 i = 0; i < c->n_cigar; ++i) {
        CigarUnit* cigar_unit = linear_alignment->add_cigar();
        cigar_unit->set_operation(kHtslibCigarToProto[bam_cigar_op(cigar[i])]);
        cigar_unit->set_operation_length(bam_cigar_oplen(cigar[i]));
      }
    }

    if (c->tid >= 0) {
      // tid >= 0 implies that the read is mapped and so has position info.
      Position* position = linear_alignment->mutable_position();
      position->set_reference_name(h->target_name[c->tid]);
      position->set_position(c->pos);
      position->set_reverse_strand(bam_is_rev(b));
    }
  }

  // Set the mates map position if the mate is not unmapped.
  // Note: According to https://samtools.github.io/hts-specs/SAMv1.pdf, RNEXT
  // field is set as '*' when the information is unavailable. htslib will
  // populate c->mtid with -1 if '*' is detected. Treat the mate as unmapped
  // even though the c->flag says otherwise.
  if (paired && !(c->flag & BAM_FMUNMAP) && c->mtid >= 0) {
    Position* mate_position = read_message->mutable_next_mate_position();

    mate_position->set_reference_name(h->target_name[c->mtid]);
    mate_position->set_position(c->mpos);
    mate_position->set_reverse_strand(bam_is_mrev(b));
  }

  // Parse out our read aux fields.
  ::nucleus::Status status = ParseAuxFields(b, options, read_message);
  if (!status.ok()) {
    // Not thread safe.
    static int counter = 0;
    if (counter++ < 1) {
      LOG(WARNING) << "Aux field parsing failure in read " << bam_get_qname(b)
                   << ": " << status;
    }
  }

  // aligned_quality may be read from aux field "OQ", therefore
  // AssignAlignedQuality function should be called after ParseAuxFields.
  status = AssignAlignedQuality(b, options, read_message);
  if (!status.ok()) {
    LOG(WARNING) << "Could not read base quality scores " << bam_get_qname(b)
                 << ": " << status;
  }

  return ::nucleus::Status();
}

// Base class for SamFullFileIterable and SamQueryIterable.
// This class implements common functionality.
class SamIterableBase : public SamIterable {
 protected:
  virtual int next_sam_record() = 0;

 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Read* out) override;

  // Base class constructor. Intializes common attrubutes.
  SamIterableBase(const SamReader* reader, htsFile* fp, bam_hdr_t* header);
  ~SamIterableBase() override;

 protected:
  htsFile* fp_;
  bam_hdr_t* header_;
  bam1_t* bam1_;
};

// Iterable class for traversing all BAM records in the file.
class SamFullFileIterable : public SamIterableBase {
 protected:
  virtual int next_sam_record();

 public:
  // Constructor is invoked via SamReader::Iterate.
  SamFullFileIterable(const SamReader* reader, htsFile* fp, bam_hdr_t* header);
};

// Iterable class for traversing BAM records returned in a query window.
class SamQueryIterable : public SamIterableBase {
 protected:
  virtual int next_sam_record();

 public:
  // Constructor will be invoked via SamReader::Query.
  SamQueryIterable(const SamReader* reader, htsFile* fp, bam_hdr_t* header,
                   hts_itr_t* iter);

  ~SamQueryIterable() override;

 private:
  hts_itr_t* iter_;
};

SamReader::SamReader(const string& reads_path, const SamReaderOptions& options,
                     htsFile* fp, bam_hdr_t* header, hts_idx_t* idx)
    : options_(options),
      fp_(fp),
      header_(header),
      idx_(idx),
      sampler_(options.downsample_fraction(), options.random_seed()) {
  CHECK(fp != nullptr) << "pointer to SAM/BAM cannot be null";
  CHECK(header_ != nullptr) << "pointer to header cannot be null";
  CHECK(options.aux_field_handling() ||
        !options.use_original_base_quality_scores())
      << "aux_field_handling must be true if use_original_quality_scores is "
         "set to true";
  bool oq_is_included = (options.aux_fields_to_keep().empty() ||
                         std::find(options.aux_fields_to_keep().begin(),
                                   options.aux_fields_to_keep().end(),
                                   kOQ) != options.aux_fields_to_keep().end());
  CHECK(oq_is_included || !options.use_original_base_quality_scores())
      << "aux_fields_to_keep must contain OQ or be empty (which means "
      << "including everything) if use_original_quality_scores is set to true.";

  const std::vector<string> header_lines_split =
      absl::StrSplit(header_->text, '\n');

  for (const string& header_line : header_lines_split) {
    const string& header_tag = header_line.substr(0, 3);
    if (header_tag == kSamHeaderTag) {
      AddHeaderLineToHeader(header_line, sam_header_);
    } else if (header_tag == kSamReferenceSequenceTag) {
      // We parse contigs separately below since they are structured by SAM
      // already.
    } else if (header_tag == kSamReadGroupTag) {
      AddReadGroupToHeader(header_line,
                           sam_header_.mutable_read_groups()->Add());
    } else if (header_tag == kSamProgramTag) {
      AddProgramToHeader(header_line, sam_header_.mutable_programs()->Add());
    } else if (header_tag == kSamCommentTag) {
      // Start at pos 4 to exclude the tab character after the tag.
      sam_header_.add_comments(header_line.substr(4));
    } else if (header_tag.empty()) {
      // Allow blank header lines
    } else {
      LOG(WARNING) << "Unrecognized SAM header type, ignoring: " << header_line;
    }
  }
  // Fill in the contig info for each contig in the sam header. Directly
  // accesses the low-level C struct because there are no indirection
  // macros/functions by htslib API.
  for (int i = 0; i < header_->n_targets; ++i) {
    nucleus::genomics::v1::ContigInfo* contig =
        sam_header_.mutable_contigs()->Add();
    contig->set_name(header_->target_name[i]);
    contig->set_n_bases(header_->target_len[i]);
    contig->set_pos_in_fasta(i);
  }
}

StatusOr<std::unique_ptr<SamReader>> SamReader::FromFile(
    const string& reads_path, const string& ref_path,
    const SamReaderOptions& options) {
  // Validate that we support the requested read requirements.
  if (options.has_read_requirements() &&
      options.read_requirements().min_base_quality_mode() !=
          nucleus::genomics::v1::ReadRequirements::UNSPECIFIED &&
      options.read_requirements().min_base_quality_mode() !=
          nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT) {
    return ::nucleus::InvalidArgument(
        absl::StrCat("Unsupported min_base_quality mode in options ",
                     options.ShortDebugString()));
  }

  htsFile* fp = hts_open_x(reads_path, "r");
  if (!fp) {
    return ::nucleus::NotFound(absl::StrCat("Could not open ", reads_path));
  }

  if (options.hts_block_size() > 0) {
    LOG(INFO) << "Setting HTS_OPT_BLOCK_SIZE to " << options.hts_block_size();
    if (hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, options.hts_block_size()) != 0)
      return ::nucleus::Unknown("Failed to set HTS_OPT_BLOCK_SIZE");
  }

  bam_hdr_t* header = sam_hdr_read(fp);
  if (header == nullptr) {
    string errmsg = absl::StrCat("bad SAM header: ", fp->fn);
    int retval = hts_close(fp);
    fp = nullptr;
    if (retval < 0) {
      return ::nucleus::Internal(
          absl::StrCat("hts_close() failed on file with ", errmsg));
    }
    return ::nucleus::Unknown(
        absl::StrCat("Could not parse file with ", errmsg));
  }

  hts_idx_t* idx = nullptr;
  if (FileTypeIsIndexable(fp->format)) {
    // TODO: use hts_idx_load after htslib upgrade.
    // This call may return null, which we will look for at Query time.
    idx = sam_index_load(fp, fp->fn);
  }

  // If we are decoding a CRAM file and the user wants to override the path to
  // the reference FASTA used to decode the CRAM, set the CRAM_OPT_REFERENCE
  // in htslib.
  if (fp->format.format == cram) {
    if (!ref_path.empty()) {
      LOG(INFO) << "Setting CRAM reference path to '" << ref_path << "'";
      if (cram_set_option(fp->fp.cram, CRAM_OPT_REFERENCE, ref_path.c_str())) {
        return ::nucleus::Unknown(absl::StrCat(
            "Failed to set the CRAM_OPT_REFERENCE value to ", ref_path));
      }
    } else {
      // If |ref_path| is empty, assumes that the reference sequence is embedded
      // in the file.
      cram_set_option(fp->fp.cram, CRAM_OPT_NO_REF, 1);
    }
  }

  return std::unique_ptr<SamReader>(
      new SamReader(reads_path, options, fp, header, idx));
}

SamReader::~SamReader() {
  if (fp_) {
    // We cannot return a value from the destructor, so the best we can do is
    // CHECK-fail if the Close() wasn't successful.
    NUCLEUS_CHECK_OK(Close());
  }
}

// Returns true if read should be returned to the client, or false otherwise.
bool SamReader::KeepRead(const nucleus::genomics::v1::Read& read) const {
  return (!options_.has_read_requirements() ||
          sam_reader_internal::ReadSatisfiesRequirements(
              read, options_.read_requirements())) &&
         // Downsample if the downsampling fraction is set.
         // Note that this can in be moved into the lower-level reader loops for
         // a slight efficiency gain (don't have to convert from bam_t to Read
         // proto but the logic to do so is much more complex than just eating
         // that cost and putting the sampling code here where it naturally fits
         // and is shared across all iteration methods.
         (options_.downsample_fraction() == 0.0 || sampler_.Keep());
}

StatusOr<std::shared_ptr<SamIterable>> SamReader::Iterate() const {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("Cannot Iterate a closed SamReader.");
  return StatusOr<std::shared_ptr<SamIterable>>(
      MakeIterable<SamFullFileIterable>(this, fp_, header_));
}

StatusOr<std::shared_ptr<SamIterable>> SamReader::Query(
    const Range& region) const {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("Cannot Query a closed SamReader.");
  if (!HasIndex()) {
    return ::nucleus::FailedPrecondition("Cannot query without an index");
  }

  const int tid = bam_name2id(header_, region.reference_name().c_str());
  if (tid < 0) {
    return ::nucleus::NotFound(
        absl::StrCat("Unknown reference_name ", region.ShortDebugString()));
  }

  // Note that query is 0-based inclusive on start and exclusive on end,
  // matching exactly the logic of our Range.
  hts_itr_t* iter = sam_itr_queryi(idx_, tid, region.start(), region.end());
  if (iter == nullptr) {
    // The region isn't valid according to sam_itr_query(), blow up.
    return ::nucleus::NotFound(
        absl::StrCat("region '", region.ShortDebugString(),
                     "' specifies an unknown reference interval"));
  }

  return StatusOr<std::shared_ptr<SamIterable>>(
      MakeIterable<SamQueryIterable>(this, fp_, header_, iter));
}

::nucleus::Status SamReader::Close() {
  if (HasIndex()) {
    hts_idx_destroy(idx_);
    idx_ = nullptr;
  }
  bam_hdr_destroy(header_);
  header_ = nullptr;
  int retval = hts_close(fp_);
  fp_ = nullptr;
  if (retval < 0) {
    return ::nucleus::Internal("hts_close() failed");
  } else {
    return ::nucleus::Status();
  }
}

// Iterable class definitions.

StatusOr<bool> SamIterableBase::Next(Read* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  // Keep reading until "reader_->KeepRead(.)"
  const SamReader* sam_reader = static_cast<const SamReader*>(reader_);
  do {
    int code = next_sam_record();
    if (code == -1) {
      return false;
    } else if (code < -1) {
      return ::nucleus::DataLoss("Failed to parse SAM record");
    }
    // Convert to proto.
    ::nucleus::Status status =
        ConvertToPb(header_, bam1_, sam_reader->options(), out);
    if (status.code() == absl::StatusCode::kAborted) {
      // "ABORT" from ConvertToPb means requirements were not met.
      continue;
    }
    NUCLEUS_RETURN_IF_ERROR(status);
  } while (!sam_reader->KeepRead(*out));
  return true;
}

SamIterableBase::SamIterableBase(const SamReader* reader, htsFile* fp,
                                 bam_hdr_t* header)
    : Iterable(reader), fp_(fp), header_(header), bam1_(bam_init1()) {}

SamIterableBase::~SamIterableBase() { bam_destroy1(bam1_); }

int SamFullFileIterable::next_sam_record() {
  // sam_read1 docs say: >= 0 on successfully reading a new record,
  // -1 on end of stream, < -1 on error.
  // Get next from file; return false if no more records to be had.
  return sam_read1(fp_, header_, bam1_);
}

SamFullFileIterable::SamFullFileIterable(const SamReader* reader, htsFile* fp,
                                         bam_hdr_t* header)
    : SamIterableBase(reader, fp, header) {}

int SamQueryIterable::next_sam_record() {
  return sam_itr_next(fp_, iter_, bam1_);
}

SamQueryIterable::~SamQueryIterable() { hts_itr_destroy(iter_); }

SamQueryIterable::SamQueryIterable(const SamReader* reader, htsFile* fp,
                                   bam_hdr_t* header, hts_itr_t* iter)
    : SamIterableBase(reader, fp, header), iter_(iter) {}

}  // namespace nucleus
