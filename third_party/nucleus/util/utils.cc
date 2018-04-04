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

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "third_party/nucleus/util/utils.h"

#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Position;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::ReadRequirements;
using nucleus::genomics::v1::Variant;

using tensorflow::strings::StrCat;
using tensorflow::StringPiece;

namespace {

// Returns a StringPiece containing the canonical base characters corresponding
// to the requested CanonicalBases in canon.
StringPiece GetCanonicalBases(const CanonicalBases canon) {
  switch (canon) {
    case CanonicalBases::ACGT:
      return "ACGT";
    case CanonicalBases::ACGTN:
      return "ACGTN";
      // We are not adding a default clause here, to explicitly make clang
      // detect the missing codes. This conversion method must stay in sync with
      // CanonicalBases enum values.
  }

  LOG(FATAL) << "Invalid CanonicalBases value" << static_cast<int>(canon);
  return "";
}

// Looks up the `variant` name in the mapping `contig_name_to_pos_in_fasta` and
// returns the corresponding pos_in_fasta. This functions assumes the variant
// name exists in the mapping.
int PosInFasta(const std::map<string, int>& contig_name_to_pos_in_fasta,
               const Variant& variant) {
  std::map<string, int>::const_iterator pos_in_fasta =
      contig_name_to_pos_in_fasta.find(variant.reference_name());
  QCHECK(pos_in_fasta != contig_name_to_pos_in_fasta.end())
      << "Reference name " << variant.reference_name()
      << " not in contig info.";
  return pos_in_fasta->second;
}

}  // namespace

bool IsCanonicalBase(const char base, const CanonicalBases canon) {
  for (const char canonical_base : GetCanonicalBases(canon)) {
    if (base == canonical_base) return true;
  }
  return false;
}

size_t FindNonCanonicalBase(StringPiece bases, const CanonicalBases canon) {
  for (size_t i = 0; i < bases.size(); i++) {
    // Meh.  This should be a static lookup table.
    if (!IsCanonicalBase(bases[i], canon)) {
      return i;
    }
  }
  return string::npos;
}

bool AreCanonicalBases(StringPiece bases, const CanonicalBases canon,
                       size_t* bad_position) {
  CHECK(!bases.empty()) << "bases cannot be empty";
  const size_t bad_pos = FindNonCanonicalBase(bases, canon);
  if (bad_pos == string::npos) return true;
  if (bad_position != nullptr) *bad_position = bad_pos;
  return false;
}

Position MakePosition(StringPiece chr, const int64 pos,
                      const bool reverse_strand) {
  Position position;
  // hm.  set_reference_name doesn't take a tensorflow::StringPiece.
  position.set_reference_name(chr.ToString());
  position.set_position(pos);
  position.set_reverse_strand(reverse_strand);
  return position;
}

Position MakePosition(const Variant& variant) {
  return MakePosition(variant.reference_name(), variant.start());
}

Range MakeRange(StringPiece chr, const int64 start, const int64 end) {
  Range range;
  range.set_reference_name(chr.ToString());
  range.set_start(start);
  range.set_end(end);
  return range;
}

Range MakeRange(const Variant& variant) {
  return MakeRange(variant.reference_name(), variant.start(), variant.end());
}

Range MakeRange(const Read& read) {
  return MakeRange(AlignedContig(read), ReadStart(read), ReadEnd(read));
}

bool RangeContains(const Range& haystack, const Range& needle) {
  return (needle.reference_name() == haystack.reference_name() &&
          needle.start() >= haystack.start() &&
          needle.end() <= haystack.end());
}

// Creates an interval string from its arguments, like chr:start-end
string MakeIntervalStr(StringPiece chr, const int64 start, const int64 end,
                       bool base_zero) {
  int offset = base_zero ? 1 : 0;
  if (start == end) {
    return StrCat(chr, ":", start + offset);
  } else {
    return StrCat(chr, ":", start + offset, "-", end + offset);
  }
}

string MakeIntervalStr(const Position& position) {
  return MakeIntervalStr(position.reference_name(), position.position(),
                         position.position(), true);
}

string MakeIntervalStr(const Range& interval) {
  return MakeIntervalStr(interval.reference_name(), interval.start(),
                         interval.end(), true);
}

string AlignedContig(const Read& read) {
  return read.has_alignment() ? read.alignment().position().reference_name()
                              : "";
}

int64 ReadStart(const Read& read) {
  return read.alignment().position().position();
}

int64 ReadEnd(const Read& read) {
  int64 position = ReadStart(read);
  for (const auto& cigar : read.alignment().cigar()) {
    switch (cigar.operation()) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::DELETE:
      case CigarUnit::SKIP:
      case CigarUnit::SEQUENCE_MISMATCH:
        position += cigar.operation_length();
        break;
      default:
        // None of the other operations change the alignment offset
        break;
    }
  }

  return position;
}

int ComparePositions(const Position& pos1, const Position& pos2) {
  int result = pos1.reference_name().compare(pos2.reference_name());
  if (result == 0) {
    result = static_cast<int>(pos1.position() - pos2.position());
  }
  return result;
}

// redacted
int ComparePositions(const Variant& variant1, const Variant& variant2) {
  return ComparePositions(MakePosition(variant1), MakePosition(variant2));
}

// True if:
// -- read is not part of a pair
// -- read is explicitly marked as properly placed by the aligner
// -- read has an unmapped mate (we only can see the next mate)
// -- read is unmapped itself
// -- read and mate are mapped to the same contig
bool IsReadProperlyPlaced(const Read& read) {
  return (read.number_reads() < 2 || read.proper_placement() ||
          read.next_mate_position().reference_name().empty() ||
          !read.has_alignment() ||
          AlignedContig(read) == read.next_mate_position().reference_name());
}

bool ReadSatisfiesRequirements(const Read& read,
                               const ReadRequirements& requirements) {
  return (requirements.keep_duplicates() || !read.duplicate_fragment()) &&
         (requirements.keep_failed_vendor_quality_checks() ||
          !read.failed_vendor_quality_checks()) &&
         (requirements.keep_secondary_alignments() ||
          !read.secondary_alignment()) &&
         (requirements.keep_supplementary_alignments() ||
          !read.supplementary_alignment()) &&
         (requirements.keep_unaligned() || read.has_alignment()) &&
         (requirements.keep_improperly_placed() ||
          IsReadProperlyPlaced(read)) &&
         (!read.has_alignment() || read.alignment().mapping_quality() >=
                                       requirements.min_mapping_quality());
}

// redacted
inline StringPiece ClippedSubstr(StringPiece s, size_t pos, size_t n) {
  pos = std::min(pos, static_cast<size_t>(s.size()));
  return s.substr(pos, n);
}

StringPiece Unquote(StringPiece input) {
  if (input.size() < 2) return input;
  char firstChar = input[0];
  char lastChar = input[input.size() - 1];
  if ((firstChar == '"' || firstChar == '\'') && (firstChar == lastChar)) {
    return ClippedSubstr(input, 1, input.size() - 2);
  } else {
    return input;
  }
}

std::map<string, int> MapContigNameToPosInFasta(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs) {
  // Create the mapping from from contig to pos_in_fasta.
  std::map<string, int> contig_name_to_pos_in_fasta;
  for (const nucleus::genomics::v1::ContigInfo& contig : contigs) {
    contig_name_to_pos_in_fasta[contig.name()] = contig.pos_in_fasta();
  }
  return contig_name_to_pos_in_fasta;
}

bool CompareVariants(const Variant& a, const Variant& b,
                     const std::map<string, int>& contig_name_to_pos_in_fasta) {
  const int pos_in_fasta_a = PosInFasta(contig_name_to_pos_in_fasta, a);
  const int pos_in_fasta_b = PosInFasta(contig_name_to_pos_in_fasta, b);
  if (pos_in_fasta_a != pos_in_fasta_b) {
    return pos_in_fasta_a < pos_in_fasta_b;
  }
  const int64 start_a = a.start();
  const int64 start_b = b.start();
  if (start_a != start_b) {
    return start_a < start_b;
  }
  return a.end() < b.end();
}

bool EndsWith(const string& s, const string& t) {
  if (t.size() > s.size()) return false;
  return std::equal(t.rbegin(), t.rend(), s.rbegin());
}

}  // namespace nucleus
