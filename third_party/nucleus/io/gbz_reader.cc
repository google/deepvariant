/*
 * Copyright 2018 Google LLC.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *o
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Implementation of gbz_reader.h
#include "third_party/nucleus/io/gbz_reader.h"

#include <algorithm>
#include <cstddef>
#include <fstream>  // IWYU pragma: keep
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "include/gbwt/metadata.h"
#include "include/gbwt/support.h"
#include "include/gbwt/utils.h"
#include "include/gbwtgraph/subgraph.h"
#include "src/include/handlegraph/types.hpp"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;

GbzReader::GbzReader(const std::string& gbz_path,
                     const std::string& sample_name)
    : sample_name_(sample_name) {
  double start = gbwt::readTimer();

  // Open GBZ file in read mode.
  std::ifstream in(gbz_path);

  // Create an empty GBZ object.
  this->gbz_ = gbwtgraph::GBZ();
  // Load the GBZ file into the GBZ object.
  gbz_.simple_sds_load(in);

  double end = gbwt::readTimer();
  LOG(INFO) << "Loading GBZ file took " << end - start << " seconds";
}

nucleus::StatusOr<std::vector<nucleus::genomics::v1::Read>> GbzReader::Query(
    const Range& range) {
  double start = gbwt::readTimer();

    const std::string& contig_name = range.reference_name();
    const size_t start_pos = range.start();
    const size_t end_pos = range.end();

    if (start_pos >= cache_start_pos_ + 300 &&
        end_pos <= std::max(cache_end_pos_ - 300, 0)) {
      // Chahe hit
      std::cerr << "[********* gbz_reader.cc line 60] Cache hit in c++:  "
              << gbwt::readTimer() - start << " seconds.\n";
      return reads_cache_;
    }


    const gbwt::Metadata& metadata = gbz_.index.metadata;

    // Get the path ids for the sample and contig
    std::vector<gbwt::size_type> path_ids = metadata.findPaths(
        metadata.sample(sample_name_), metadata.contig(contig_name));

    if (path_ids.empty()) {
      return nucleus::NotFound(absl::StrCat(
          "Pangenome path ids not found for pangenome sample name: ",
          sample_name_));
    }
    handlegraph::path_handle_t path =
        gbz_.graph.path_to_handle(path_ids.front());

    gbwtgraph::SubgraphQuery query = gbwtgraph::SubgraphQuery::path_interval(
        path, start_pos, end_pos, /*context=*/ 1000,
        gbwtgraph::SubgraphQuery::HaplotypeOutput::all_haplotypes);

    auto path_index = std::make_unique<gbwtgraph::PathIndex>(gbz_);

    gbwtgraph::Subgraph subgraph(gbz_, path_index.get(), query);

    const std::vector<nucleus::genomics::v1::Read>& reads =
        GetReadsFromSubgraph(subgraph, gbz_);

    std::cerr << "Query in c++:  "
              << gbwt::readTimer() - start << " seconds.\n";

    updateCache(reads);

    return nucleus::StatusOr<std::vector<nucleus::genomics::v1::Read>>(reads);
}

nucleus::StatusOr<std::shared_ptr<SamIterable>> GbzReader::Iterate() const {
  return nucleus::Unimplemented("Iterate is not implemented for GbzReader.");
}

void GbzReader::updateCache(
    const std::vector<nucleus::genomics::v1::Read>& reads) {
  reads_cache_ = reads;
  if (reads.empty()) {
    cache_start_pos_ = 0;
    cache_end_pos_ = 0;
    return;
  }

  cache_start_pos_ = std::max_element(
                         reads_cache_.begin(), reads_cache_.end(),
                         [](const nucleus::genomics::v1::Read& a,
                            const nucleus::genomics::v1::Read& b) {
                           return a.alignment().position().position() <
                                  b.alignment().position().position();
                         })
                         ->alignment()
                         .position()
                         .position();

  auto last_read =
      std::min_element(reads_cache_.begin(), reads_cache_.end(),
                       [](const nucleus::genomics::v1::Read& a,
                          const nucleus::genomics::v1::Read& b) {
                         return a.alignment().position().position() +
                                    a.aligned_sequence().size() >
                                b.alignment().position().position() +
                                    b.aligned_sequence().size();
                       });
  cache_end_pos_ = last_read->alignment().position().position() +
                   last_read->aligned_sequence().size();
}

std::vector<nucleus::genomics::v1::Read> GbzReader::GetReadsFromSubgraph(
    const gbwtgraph::Subgraph& subgraph, const gbwtgraph::GBZ& gbz) {
  // W-lines: reference path
  std::string contig_name = "unknown";
  int start = 0;
  if (subgraph.reference_path < subgraph.paths.size()) {
    // const gbwt::vector_type& path =
    // subgraph.paths[subgraph.reference_path];
    gbwt::FullPathName path_name = subgraph.reference_path_name(gbz);
    contig_name = path_name.contig_name;
    start = path_name.offset;
  }
  std::vector<nucleus::genomics::v1::Read> reads;

  // W-lines: anonymous haplotypes
  size_t haplotype_id = 1;
  for (size_t i = 0; i < subgraph.paths.size(); i++) {
    if (i == subgraph.reference_path) {
      continue;
    }
    const gbwt::vector_type& path = subgraph.paths[i];
    const std::string* cigar = subgraph.cigar(i);
    // get a vector of cigar elements
    const std::vector<std::string> cigar_elements =
        GbzReader::GetCigarElements(*cigar);
    // get the sequence of the haplotype
    const std::string bases = GbzReader::GetBases(path, gbz);
    // make a Read object for the haplotype and its alignment
    nucleus::genomics::v1::Read read =
        MakeRead(contig_name, start, bases, cigar_elements,
                 std::string("haplotype_") + std::to_string(haplotype_id));
    // add Read object to the vector
    reads.push_back(read);
    haplotype_id++;
  }
  return reads;
}

nucleus::genomics::v1::Read GbzReader::MakeRead(
    const std::string& chr, const int start, const string& bases,
    const std::vector<string>& cigar_elements, const string& read_name) {
  Read read;
  read.set_fragment_name(read_name);
  read.set_aligned_sequence(bases);
  read.set_number_reads(2);
  read.set_proper_placement(true);

  for (size_t i = 0; i < bases.length(); ++i) {
    read.add_aligned_quality(30);
  }

  genomics::v1::LinearAlignment& aln = *read.mutable_alignment();
  aln.set_mapping_quality(90);
  *aln.mutable_position() = nucleus::MakePosition(chr, start);
  for (const auto& cigar_str : cigar_elements) {
    genomics::v1::CigarUnit& cigar = *aln.add_cigar();
    cigar.set_operation(ParseCigarOpStr(cigar_str.back()));
    const string& cigar_len_str = cigar_str.substr(0, cigar_str.length() - 1);
    int value = atoi(cigar_len_str.c_str());  // NOLINT: safe here.
    cigar.set_operation_length(value);
  }

  return read;
}

genomics::v1::CigarUnit_Operation GbzReader::ParseCigarOpStr(const char op) {
  switch (op) {
    case 'M':
      return genomics::v1::CigarUnit::ALIGNMENT_MATCH;
    case '=':
      return genomics::v1::CigarUnit::SEQUENCE_MATCH;
    case 'X':
      return genomics::v1::CigarUnit::SEQUENCE_MISMATCH;
    case 'I':
      return genomics::v1::CigarUnit::INSERT;
    case 'D':
      return genomics::v1::CigarUnit::DELETE;
    case 'S':
      return genomics::v1::CigarUnit::CLIP_SOFT;
    case 'P':
      return genomics::v1::CigarUnit::PAD;
    case 'H':
      return genomics::v1::CigarUnit::CLIP_HARD;
    case 'N':
      return genomics::v1::CigarUnit::SKIP;
    default:
      LOG(FATAL) << "Unexpected cigar op " << op;
  }
}
std::vector<std::string> GbzReader::GetCigarElements(const std::string& input) {
  const std::regex pattern(R"(\d+[MID=XHS])");
  // Vector to store the matched results
  std::vector<std::string> results;

  // Using a regex_iterator to find all matches
  auto matches_begin =
      std::sregex_iterator(input.begin(), input.end(), pattern);
  auto matches_end = std::sregex_iterator();

  // Iterate over all matches and store them in the results vector
  for (std::sregex_iterator i = matches_begin; i != matches_end; ++i) {
    const std::smatch& match = *i;
    results.push_back(match.str());
  }
  return results;
}

std::string GbzReader::GetBases(const gbwt::vector_type& path,
                                const gbwtgraph::GBZ& gbz) {
  std::string bases("");
  for (gbwt::node_type node : path) {
    handlegraph::handle_t handle = gbz.graph.get_handle(gbwt::Node::id(node));
    if (gbwt::Node::is_reverse(node)) {
      bases += GbzReader::GetReverseComplement(gbz.graph.get_sequence(handle));
    } else {
      bases += gbz.graph.get_sequence(handle);
    }
  }
  return bases;
}

std::string GbzReader::GetReverseComplement(const std::string& sequence) {
  std::string complement;
  // Iterate through each character in the sequence
  for (char nucleotide : sequence) {
    // Convert the character to its complement
    switch (nucleotide) {
      case 'A':
        complement += 'T';
        break;
      case 'T':
        complement += 'A';
        break;
      case 'C':
        complement += 'G';
        break;
      case 'G':
        complement += 'C';
        break;
      default:
        // If the character is not a valid nucleotide, add it unchanged
        complement += nucleotide;
        break;
    }
  }
  // Reverse the complement string
  std::reverse(complement.begin(), complement.end());

  return complement;
}

}  // namespace nucleus
