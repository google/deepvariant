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

#include "third_party/nucleus/testing/test_utils.h"

#include <vector>

#include "absl/strings/str_join.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

using nucleus::MakePosition;
using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::CigarUnit_Operation;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Read;

namespace {

// Gets the CigarUnit_Operation corresponding to the single character
// representation of that CIGAR option.  For example, 'M' returns
// CigarUnit::ALIGNMENT_MATCH.
CigarUnit_Operation parse_cigar_op_str(const char op) {
  switch (op) {
    case 'M':
      return CigarUnit::ALIGNMENT_MATCH;
    case '=':
      return CigarUnit::SEQUENCE_MATCH;
    case 'X':
      return CigarUnit::SEQUENCE_MISMATCH;
    case 'I':
      return CigarUnit::INSERT;
    case 'D':
      return CigarUnit::DELETE;
    case 'S':
      return CigarUnit::CLIP_SOFT;
    case 'P':
      return CigarUnit::PAD;
    case 'H':
      return CigarUnit::CLIP_HARD;
    case 'N':
      return CigarUnit::SKIP;
    default:
      LOG(FATAL) << "Unexpected cigar op " << op;
  }
}

// TODO: only supports UNIX paths.
template<class T>
string JoinPaths(std::initializer_list<T> paths) {
  return absl::StrJoin(paths, "/");
}

}  // namespace

// Simple getter for test files in the right testdata path.
string GetTestData(absl::string_view path, absl::string_view test_data_dir) {
  // const string test_srcdir = getenv("TEST_SRCDIR");
  // const char* test_workspace = getenv("TEST_WORKSPACE");
  absl::string_view test_srcdir = getenv("TEST_SRCDIR");
  absl::string_view test_workspace = getenv("TEST_WORKSPACE");
  test_workspace = !test_workspace.empty() ? test_workspace : kDefaultWorkspace;
  return JoinPaths({test_srcdir, test_workspace, test_data_dir, path});
}

string MakeTempFile(absl::string_view filename) {
  absl::string_view test_tmpdir = getenv("TEST_TMPDIR");
  return JoinPaths({test_tmpdir, filename});
}

std::vector<nucleus::genomics::v1::ContigInfo> CreateContigInfos(
    const std::vector<string>& names, const std::vector<int>& positions) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs;
  QCHECK_EQ(names.size(), positions.size());
  nucleus::genomics::v1::ContigInfo contig;
  for (size_t i = 0; i < names.size(); i++) {
    contig.set_name(names[i]);
    contig.set_pos_in_fasta(positions[i]);
    contigs.push_back(contig);
  }
  return contigs;
}

Read MakeRead(const string& chr,
              const int start, const string& bases,
              const std::vector<string>& cigar_elements,
              const string& read_name) {
  Read read;
  read.set_fragment_name(read_name);
  read.set_aligned_sequence(bases);
  read.set_number_reads(2);
  read.set_proper_placement(true);

  for (size_t i = 0; i < bases.length(); ++i) {
    read.add_aligned_quality(30);
  }

  LinearAlignment& aln = *read.mutable_alignment();
  aln.set_mapping_quality(90);
  *aln.mutable_position() = MakePosition(chr, start);
  for (const auto& cigar_str : cigar_elements) {
    CigarUnit& cigar = *aln.add_cigar();
    cigar.set_operation(parse_cigar_op_str(cigar_str.back()));
    const string& cigar_len_str = cigar_str.substr(0, cigar_str.length() - 1);
    int value = atoi(cigar_len_str.c_str());  // NOLINT: safe here.
    cigar.set_operation_length(value);
  }

  return read;
}

std::vector<CigarUnit> MakeCigar(
    const std::vector<std::string>& cigar_elements) {
  std::vector<CigarUnit> cigar;
  for (const auto& cigar_str : cigar_elements) {
    CigarUnit cigar_unit = CigarUnit();
    cigar_unit.set_operation(parse_cigar_op_str(cigar_str.back()));
    const string& cigar_len_str = cigar_str.substr(0, cigar_str.length() - 1);
    int value = atoi(cigar_len_str.c_str());  // NOLINT: safe here.
    cigar_unit.set_operation_length(value);
    cigar.push_back(cigar_unit);
  }
  return cigar;
}

bool IsGzipped(absl::string_view file_contents) {
  const char gzip_magic[2] = {'\x1f', '\x8b'};
  return (file_contents.size() >= 2 && file_contents[0] == gzip_magic[0] &&
          file_contents[1] == gzip_magic[1]);
}



}  // namespace nucleus
