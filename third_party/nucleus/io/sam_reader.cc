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

// Implementation of sam_reader.h
#include "third_party/nucleus/io/sam_reader.h"

#include "google/protobuf/repeated_field.h"
#include "htslib/hts.h"
#include "htslib/hts_endian.h"
#include "htslib/sam.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/index.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/strings/numbers.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

namespace tf = tensorflow;

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::CigarUnit_Operation;
using nucleus::genomics::v1::Position;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::SamHeader;
using nucleus::genomics::v1::SamReaderOptions;
using std::vector;

using tensorflow::WARNING;
using tensorflow::int32;

using google::protobuf::RepeatedField;
using tensorflow::strings::StrCat;

namespace {
constexpr char kSamHeaderTag[] = "@HD";
constexpr char kSamReferenceSequenceTag[] = "@SQ";
constexpr char kSamReadGroupTag[] = "@RG";
constexpr char kSamProgramTag[] = "@PG";
constexpr char kSamCommentTag[] = "@CO";

void AddHeaderLineToHeader(const string& line, SamHeader& header) {
  static constexpr char kVersionTag[] = "VN:";
  static constexpr char kSortingOrderTag[] = "SO:";
  static constexpr char kAlignmentGroupingTag[] = "GO:";
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

  for (const string& token : tensorflow::str_util::Split(line, '\t')) {
    if (token == kSamHeaderTag) continue;
    const string tag = token.substr(0, tagLen);
    const string value = token.substr(tagLen);
    if (tag == kVersionTag) {
      header.set_format_version(value);
    } else if (tag == kSortingOrderTag) {
      const auto& it = sorting_order_map.find(value);
      if (it == sorting_order_map.end()) {
        LOG(WARNING) << "Unknown sorting order, defaulting to unknown: "
                     << line;
        header.set_sorting_order(SamHeader::UNKNOWN);
      } else {
        header.set_sorting_order(it->second);
      }
    } else if (tag == kAlignmentGroupingTag) {
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
  for (const string& token : tensorflow::str_util::Split(line, '\t')) {
    if (token == kSamReadGroupTag) continue;
    const string tag = token.substr(0, tagLen);
    const string value = token.substr(tagLen);
    if (tag == "ID:") {
      readgroup->set_name(value);
    } else if (tag == "CN:") {
      readgroup->set_sequencing_center(value);
    } else if (tag == "DS:") {
      readgroup->set_description(value);
    } else if (tag == "DT:") {
      readgroup->set_date(value);
    } else if (tag == "FO:") {
      readgroup->set_flow_order(value);
    } else if (tag == "KS:") {
      readgroup->set_key_sequence(value);
    } else if (tag == "LB:") {
      readgroup->set_library_id(value);
    } else if (tag == "PG:") {
      readgroup->add_program_ids(value);
    } else if (tag == "PI:") {
      int size;
      tensorflow::strings::safe_strto32(value, &size);
      readgroup->set_predicted_insert_size(size);
    } else if (tag == "PL:") {
      readgroup->set_platform(value);
    } else if (tag == "PM:") {
      readgroup->set_platform_model(value);
    } else if (tag == "PU:") {
      readgroup->set_platform_unit(value);
    } else if (tag == "SM:") {
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
  for (const string& token : tensorflow::str_util::Split(line, '\t')) {
    if (token == kSamProgramTag) continue;
    const string tag = token.substr(0, tagLen);
    const string value = token.substr(tagLen);
    if (tag == "ID:") {
      program->set_id(value);
    } else if (tag == "PN:") {
      program->set_name(value);
    } else if (tag == "CL:") {
      program->set_command_line(value);
    } else if (tag == "PP:") {
      program->set_prev_program_id(value);
    } else if (tag == "DS:") {
      program->set_description(value);
    } else if (tag == "VN:") {
      program->set_version(value);
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------
//
// Reader for SAM/BAM/CRAM etc formats containing NGS reads supported by htslib.
//
// -----------------------------------------------------------------------------

// Array mapping htslib BAM constants (in comment) to proto CigarUnit enum
// values.
const CigarUnit_Operation kHtslibCigarToProto[] = {
// #define BAM_CMATCH      0
  CigarUnit::ALIGNMENT_MATCH,
// #define BAM_CINS        1
  CigarUnit::INSERT,
// #define BAM_CDEL        2
  CigarUnit::DELETE,
// #define BAM_CREF_SKIP   3
  CigarUnit::SKIP,
// #define BAM_CSOFT_CLIP  4
  CigarUnit::CLIP_SOFT,
// #define BAM_CHARD_CLIP  5
  CigarUnit::CLIP_HARD,
// #define BAM_CPAD        6
  CigarUnit::ALIGNMENT_MATCH,
// #define BAM_CEQUAL      7
  CigarUnit::SEQUENCE_MATCH,
// #define BAM_CDIFF       8
  CigarUnit::SEQUENCE_MISMATCH,
// #define BAM_CBACK       9
  CigarUnit::OPERATION_UNSPECIFIED,
};

// Gets the size in bytes for a SAM/BAM aux tag based on their declared type.
//
// Based on code from htslib/sam.c which isn't exported from htslib. Returns
// a value < 0 if type isn't one of the expected types from htslib.
static inline int HtslibAuxSize(uint8_t type) {
  switch (type) {
    case 'A': case 'c': case 'C':
      return 1;
    case 's': case 'S':
      return 2;
    case 'f': case 'i': case 'I':
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
//  From https://samtools.github.io/hts-specs/SAMv1.pdf
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
tf::Status ParseAuxFields(const bam1_t* b, const SamReaderOptions& options,
                          Read* read_message) {
  if (options.aux_field_handling() != SamReaderOptions::PARSE_ALL_AUX_FIELDS) {
    return tf::Status::OK();
  }

  uint8_t* s = bam_get_aux(b);
  const uint8_t* end = b->data + b->l_data;
  while (end - s >= 4) {
    // Each block is encoded like (each element is a byte):
    // [tag char 1, tag char 2, type byte, ...]
    // where the ... contents depends on the 2-character tag and type.
    const string tag = string(reinterpret_cast<char*>(s), 2);
    s += 2;
    const uint8_t type = *s++;
    switch (type) {
      // An 'A' is just a single character string.
      case 'A': {
        // Safe since we know s is at least 4 bytes from the end.
        const string value = string(reinterpret_cast<char*>(s), 1);
        SetInfoField(tag, value, read_message);
        s += 1;
      } break;
      // These are all different byte-sized integers.
      case 'C': case 'c': case 'S': case 's': case 'I': case 'i': {
        const int size = HtslibAuxSize(type);
        if (size < 0 || end - s < size)
          return tf::errors::DataLoss("Malformed tag " + tag);
        errno = 0;
        const int value = bam_aux2i(s - 1);
        if (value == 0 && errno == EINVAL)
          return tf::errors::DataLoss("Malformed tag " + tag);
        SetInfoField(tag, value, read_message);
        s += size;
      } break;
      // A 4-byte floating point.
      case 'f': {
        if (end - s < 4) return tf::errors::DataLoss("Malformed tag " + tag);
        const float value = le_to_float(s);
        SetInfoField(tag, value, read_message);
        s += 4;
      } break;
      // Z and H are null-terminated strings.
      case 'Z': case 'H': {
        char* value = reinterpret_cast<char*>(s);
        for (; s < end && *s; ++s) {}  // Loop to the end.
        if (s >= end) return tf::errors::DataLoss("Malformed tag " + tag);
        s++;
        // The H hex tag is not really used and likely deprecated (see:
        // https://sourceforge.net/p/samtools/mailman/message/28274509/
        // so we are explicitly skipping them here.
        if (type == 'Z') SetInfoField(tag, value, read_message);
      } break;
      // B is an array of atomic types (strings, ints, floats).
      case 'B': {
        const uint8_t sub_type = *s++;
        const int element_size = HtslibAuxSize(sub_type);
        if (element_size < 0)
          return tf::errors::DataLoss("element_size == 0 for tag " + tag);
        // Prevents us from reading off the end of our buffer with le_to_u32.
        if (end - s < 4)
          return tf::errors::DataLoss("data too short for tag " + tag);
        const int n_elements = le_to_u32(s);
        if (n_elements == 0) return tf::errors::DataLoss("n_elements is zero");
        // redacted
        // We need to skip len * element size in bytes to clear the array
        // itself, and 4 more bytes for n_elements int that occurs before the
        // array.
        s += 4 + n_elements * element_size;
      } break;
      default: {
        return tf::errors::DataLoss("Unknown tag " + tag);
      }
    }
  }

  // Everything parsed correctly, so we return OK.
  return tf::Status::OK();
}

tf::Status ConvertToPb(const bam_hdr_t* h, const bam1_t* b,
                       const SamReaderOptions& options, Read* read_message) {
  CHECK(h != nullptr) << "BAM header cannot be null";
  CHECK(b != nullptr) << "BAM record cannot be null";
  CHECK(read_message != nullptr) << "Read record cannot be null";

  read_message->Clear();

  const bam1_core_t *c = &b->core;

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

  if (c->l_qseq) {
    // Convert the seq and qual fields if they are present.
    string* read_seq = read_message->mutable_aligned_sequence();
    read_seq->reserve(c->l_qseq);
    uint8_t* seq = bam_get_seq(b);  // seq is stored as 8-bit offsets.
    for (int i = 0; i < c->l_qseq; ++i) {
      // Convert the offsets to upper case characters by their offset in
      // the constant seq_nt16_str from htslib.
      read_seq->push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }

    // Convert the qual field.
    uint8_t* quals = bam_get_qual(b);
    if (quals[0] != 0xff) {  // Not missing
      // redacted
      RepeatedField<int32>* quality = read_message->mutable_aligned_quality();
      quality->Reserve(c->l_qseq);
      for (int i = 0; i < c->l_qseq; ++i) {
        quality->Add(quals[i]);
      }
    }
  }

  if (!(c->flag & BAM_FUNMAP)) {
    // If the read is mapped, set the mapping information in read_message.
    auto* linear_alignment = read_message->mutable_alignment();
    linear_alignment->set_mapping_quality(c->qual);

    if (c->n_cigar) {  // Convert our Cigar.
      uint32_t* cigar = bam_get_cigar(b);
      for (uint32_t i = 0; i < c->n_cigar; ++i) {
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
  if (paired && !(c->flag & BAM_FMUNMAP)) {
    Position* mate_position = read_message->mutable_next_mate_position();
    if (c->mtid < 0)
      return tf::errors::DataLoss(
          StrCat("Expected mtid >= 0 as mate is supposedly mapped: ",
                 read_message->ShortDebugString()));
    mate_position->set_reference_name(h->target_name[c->mtid]);
    mate_position->set_position(c->mpos);
    mate_position->set_reverse_strand(bam_is_mrev(b));
  }

  // Parse out our read aux fields.
  tf::Status status = ParseAuxFields(b, options, read_message);
  if (!status.ok()) {
    // Not thread safe.
    static int counter = 0;
    if (counter++ < 1) {
      LOG(WARNING) << "Aux field parsing failure in read "
                   << bam_get_qname(b) << ": " << status;
    }
  }

  return tf::Status::OK();
}

// Iterable class for traversing all BAM records in the file.
class SamFullFileIterable : public SamIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Read* out) override;

  // Constructor is invoked via SamReader::Iterate.
  SamFullFileIterable(const SamReader* reader, htsFile* fp, bam_hdr_t* header);
  ~SamFullFileIterable() override;

 private:
  htsFile* fp_;
  bam_hdr_t* header_;
  bam1_t* bam1_;
};


// Iterable class for traversing BAM records returned in a query window.
class SamQueryIterable : public SamIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Read* out) override;

  // Constructor will be invoked via SamReader::Query.
  SamQueryIterable(const SamReader* reader,
                   htsFile* fp,
                   bam_hdr_t* header,
                   hts_itr_t* iter);

  ~SamQueryIterable() override;

 private:
  htsFile* fp_;
  bam_hdr_t* header_;
  hts_itr_t* iter_;
  bam1_t* bam1_;
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

  const std::vector<string> header_lines_split =
      tensorflow::str_util::Split(header_->text, '\n');

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
      sam_header_.add_comments(header_line);
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
    const string& reads_path, const SamReaderOptions& options) {
  // Validate that we support the requested read requirements.
  if (options.has_read_requirements() &&
      options.read_requirements().min_base_quality_mode() !=
          nucleus::genomics::v1::ReadRequirements::UNSPECIFIED &&
      options.read_requirements().min_base_quality_mode() !=
          nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT) {
    return tf::errors::InvalidArgument(
        StrCat("Unsupported min_base_quality mode in options ",
               options.ShortDebugString()));
  }

  htsFile* fp = hts_open_x(reads_path.c_str(), "r");
  if (!fp) {
    return tf::errors::NotFound(StrCat("Could not open ", reads_path));
  }

  if (options.hts_block_size() > 0) {
    LOG(INFO) << "Setting HTS_OPT_BLOCK_SIZE to " << options.hts_block_size();
    if (hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, options.hts_block_size()) != 0)
      return tf::errors::Unknown(StrCat("Failed to set HTS_OPT_BLOCK_SIZE"));
  }

  bam_hdr_t* header = sam_hdr_read(fp);
  if (header == nullptr)
    return tf::errors::Unknown(StrCat("Couldn't parse header for ", fp->fn));

  hts_idx_t* idx = nullptr;
  if (options.index_mode() ==
      nucleus::genomics::v1::IndexHandlingMode::INDEX_BASED_ON_FILENAME) {
    // redacted
    idx = sam_index_load(fp, fp->fn);
    if (idx == nullptr) {
      return tf::errors::NotFound(StrCat("No index found for ", fp->fn));
    }
  }

  return std::unique_ptr<SamReader>(
      new SamReader(reads_path, options, fp, header, idx));
}

SamReader::~SamReader() {
  if (fp_) {
    // We cannot return a value from the destructor, so the best we can do is
    // CHECK-fail if the Close() wasn't successful.
    TF_CHECK_OK(Close());
  }
}

// Returns true if read should be returned to the client, or false otherwise.
bool SamReader::KeepRead(const nucleus::genomics::v1::Read& read) const {
  return (!options_.has_read_requirements() ||
          ReadSatisfiesRequirements(read, options_.read_requirements())) &&
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
    return tf::errors::FailedPrecondition("Cannot Iterate a closed SamReader.");
  return StatusOr<std::shared_ptr<SamIterable>>(
      MakeIterable<SamFullFileIterable>(this, fp_, header_));
}

StatusOr<std::shared_ptr<SamIterable>> SamReader::Query(
    const Range& region) const {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition("Cannot Query a closed SamReader.");
  if (!HasIndex()) {
    return tf::errors::FailedPrecondition("Cannot query without an index");
  }

  const int tid = bam_name2id(header_, region.reference_name().c_str());
  if (tid < 0) {
    return tf::errors::NotFound(
        StrCat("Unknown reference_name ", region.ShortDebugString()));
  }

  // Note that query is 0-based inclusive on start and exclusive on end,
  // matching exactly the logic of our Range.
  hts_itr_t* iter = sam_itr_queryi(idx_, tid, region.start(), region.end());
  if (iter == nullptr) {
    // The region isn't valid according to sam_itr_query(), blow up.
    return tf::errors::NotFound(
        StrCat("region '", region.ShortDebugString(),
               "' specifies an unknown reference interval"));
  }

  return StatusOr<std::shared_ptr<SamIterable>>(
      MakeIterable<SamQueryIterable>(this, fp_, header_, iter));
}


tf::Status SamReader::Close() {
  if (HasIndex()) {
    hts_idx_destroy(idx_);
    idx_ = nullptr;
  }
  bam_hdr_destroy(header_);
  header_ = nullptr;
  int retval = hts_close(fp_);
  fp_ = nullptr;
  if (retval < 0) {
    return tf::errors::Internal("hts_close() failed");
  } else {
    return tf::Status::OK();
  }
}


// Iterable class definitions.

StatusOr<bool> SamFullFileIterable::Next(Read* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  // Keep reading until "reader_->KeepRead(.)"
  const SamReader* sam_reader = static_cast<const SamReader*>(reader_);
  do {
    // sam_read1 docs say: >= 0 on successfully reading a new record,
    // -1 on end of stream, < -1 on error.
    // Get next from file; return false if no more records to be had.
    int code = sam_read1(fp_, header_, bam1_);
    if (code == -1) {
      return false;
    } else if (code < -1) {
      return tf::errors::DataLoss("Failed to parse SAM record");
    }
    // Convert to proto.
    TF_RETURN_IF_ERROR(ConvertToPb(header_, bam1_, sam_reader->options(), out));
  } while (!sam_reader->KeepRead(*out));
  return true;
}

SamFullFileIterable::~SamFullFileIterable() {
  bam_destroy1(bam1_);
}

SamFullFileIterable::SamFullFileIterable(const SamReader* reader,
                                         htsFile* fp,
                                         bam_hdr_t* header)
    : Iterable(reader),
      fp_(fp),
      header_(header),
      bam1_(bam_init1())
{}

// redacted
// class that only differs in sam_itr_next vs sam_read1 calls.
StatusOr<bool> SamQueryIterable::Next(Read* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  // Keep reading until "reader_->KeepRead(.)"
  const SamReader* sam_reader = static_cast<const SamReader*>(reader_);
  do {
    // Get next in query window; return false if no more records.
    int code = sam_itr_next(fp_, iter_, bam1_);
    if (code == -1) {
      return false;
    } else if (code < -1) {
      return tf::errors::DataLoss("Failed to parse SAM record");
    }
    // Convert to proto.
    TF_RETURN_IF_ERROR(ConvertToPb(header_, bam1_, sam_reader->options(), out));
  } while (!sam_reader->KeepRead(*out));
  return true;
}

SamQueryIterable::~SamQueryIterable() {
  bam_destroy1(bam1_);
  hts_itr_destroy(iter_);
}

SamQueryIterable::SamQueryIterable(const SamReader* reader,
                                   htsFile* fp,
                                   bam_hdr_t* header,
                                   hts_itr_t* iter)
    : Iterable(reader),
      fp_(fp),
      header_(header),
      iter_(iter),
      bam1_(bam_init1())
{}

}  // namespace nucleus
