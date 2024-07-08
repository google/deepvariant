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

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "third_party/nucleus/util/utils.h"

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/substitute.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace nucleus {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Position;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::Variant;

using absl::string_view;
using absl::Substitute;

namespace {

// Returns a StringPiece containing the canonical base characters corresponding
// to the requested CanonicalBases in canon.
string_view GetCanonicalBases(const CanonicalBases canon) {
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
  return GetCanonicalBases(canon).find(base) != string_view::npos;
}

size_t FindNonCanonicalBase(string_view bases, const CanonicalBases canon) {
  for (size_t i = 0; i < bases.size(); i++) {
    // Meh.  This should be a static lookup table.
    if (!IsCanonicalBase(bases[i], canon)) {
      return i;
    }
  }
  return string::npos;
}

bool AreCanonicalBases(string_view bases, const CanonicalBases canon,
                       size_t* bad_position) {
  CHECK(!bases.empty()) << "bases cannot be empty";
  const size_t bad_pos = FindNonCanonicalBase(bases, canon);
  if (bad_pos == string::npos) return true;
  if (bad_position != nullptr) *bad_position = bad_pos;
  return false;
}

Position MakePosition(string_view chr, const int64 pos,
                      const bool reverse_strand) {
  Position position;
  position.set_reference_name(string(chr));  // TODO
  position.set_position(pos);
  position.set_reverse_strand(reverse_strand);
  return position;
}

Position MakePosition(const Variant& variant) {
  return MakePosition(variant.reference_name(), variant.start());
}

Range MakeRange(string_view chr, const int64 start, const int64 end) {
  Range range;
  range.set_reference_name(string(chr));  // TODO
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

void ReadRangePython(
    const nucleus::ConstProtoPtr<const ::nucleus::genomics::v1::Read>&
        read_wrapped,
    nucleus::EmptyProtoPtr<::nucleus::genomics::v1::Range> range_wrapped) {
  const Read& read = *read_wrapped.p_;
  Range* range = range_wrapped.p_;
  range->set_reference_name(read.alignment().position().reference_name());
  range->set_start(ReadStart(read));
  range->set_end(ReadEnd(read));
}

bool RangeContains(const Range& haystack, const Range& needle) {
  return (needle.reference_name() == haystack.reference_name() &&
          needle.start() >= haystack.start() &&
          needle.end() <= haystack.end());
}

bool RangesContainVariant(
    const std::vector<nucleus::genomics::v1::Range>& ranges,
    const Variant& variant) {
  for (const auto& range : ranges) {
    if (range.reference_name() == variant.reference_name() &&
        range.start() <= variant.start() && range.end() > variant.start()) {
      return true;
    }
  }
  return false;
}

bool ReadOverlapsRegion(const ::nucleus::genomics::v1::Read& read,
                        const ::nucleus::genomics::v1::Range& range) {
  // Equivalent code in python from ranges.py:
  //
  // return (i1.reference_name == i2.reference_name and i1.end > i2.start and
  //         i1.start < i2.end)
  //
  // Here i1 is range and i2 is the range implied from the read.
  return
      // This is the cheapest calculation as read start is cheap to determine.
      range.end() > ReadStart(read) &&
      // Next we check read end, which is slightly more expensive as we need to
      // compute the end from the cigar.
      range.start() < ReadEnd(read) &&
      // Finally we compute if the reference_names are the same.
      range.reference_name() == AlignedContig(read);
}

// Creates an interval string from its arguments, like chr:start-end
string MakeIntervalStr(string_view chr, const int64 start, const int64 end,
                       bool base_zero) {
  int offset = base_zero ? 1 : 0;
  if (start == end) {
    // TODO: remove string conversion
    return Substitute("$0:$1", string(chr), start + offset);
  } else {
    // TODO: remove string conversion
    return Substitute("$0:$1-$2", string(chr), start + offset, end + offset);
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

// TODO: should compare ranges, implement compare range
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

inline string_view ClippedSubstr(string_view s, size_t pos, size_t n) {
  pos = std::min(pos, static_cast<size_t>(s.size()));
  return s.substr(pos, n);
}

string_view Unquote(string_view input) {
  if (input.size() < 2) return input;
  char firstChar = input[0];
  char lastChar = input[input.size() - 1];
  if ((firstChar == '"' || firstChar == '\'') && (firstChar == lastChar)) {
    return nucleus::ClippedSubstr(input, 1, input.size() - 2);
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
