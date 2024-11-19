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

#include "boost/interprocess/managed_shared_memory.hpp"  // NOLINT
#include "boost/interprocess/shared_memory_object.hpp"  // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp"  // NOLINT


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
#include "include/gbwt/utils.h"
#include "include/gbwtgraph/subgraph.h"
#include "include/gbwtgraph/gbz.h"
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
namespace bi = boost::interprocess;

void GbzReader::create_or_open_shared_memory(){
  if (this->shared_memory_ != nullptr){
    LOG(WARNING) << "shared memory is not null so not opening it again\n";
    return;
  }
  if (this->shared_memory_name_.empty()){
    LOG(FATAL) << "shared memory should have a non-empty name\n";
    return;
  }
  if (this->shared_memory_size_gb_ <= 0){
    LOG(FATAL) << "shared_memory should have a non-empty name\n";
    return;
  }

  if (this->create_shared_memory_) {
    bi::shared_memory_object::remove(this->shared_memory_name_.c_str());
    LOG(INFO) << "Creating shared memory: " << this->shared_memory_name_;
    // Create shared memory
    this->shared_memory_ = new bi::managed_shared_memory(
        bi::create_only,
        this->shared_memory_name_.c_str(),
        shared_memory_size_gb_ * 1e9);
    // set the process countdown to the number of processes that will use the
    // shared memory
    // lock the mutex to make sure no other process is modifying the countdown
    bi::named_mutex mtx_process_countdown(
        bi::open_or_create,
        "PROCESS_COUNTDOWN_MUTEX");
    bi::scoped_lock<bi::named_mutex> lock(mtx_process_countdown);
    int *process_countdown = this->shared_memory_->find_or_construct<int>
        ("PROCESS_COUNTDOWN")(0);
    *process_countdown = num_processes_;
    LOG(INFO) <<
        "Number of expected processes (excluding the one that created"
        " the shared memory) using "
        "shared memory is set to:\n" <<
        *process_countdown << "\n";
  } else {
    LOG(INFO) << "Opening shared memory " << this->shared_memory_name_;
    this->shared_memory_ = new bi::managed_shared_memory(
        bi::open_only,
        this->shared_memory_name_.c_str());
  }
}

void GbzReader::close_shared_memory(){
  if (this->shared_memory_ != nullptr){
    // lock the mutex to make sure no other process is modifying the countdown
    bi::named_mutex mtx_process_countdown(
        bi::open_or_create,
        "PROCESS_COUNTDOWN_MUTEX");
    bi::scoped_lock<bi::named_mutex> lock(mtx_process_countdown);
    int *process_countdown = this->shared_memory_->find_or_construct<int>
        ("PROCESS_COUNTDOWN")(0);
    if ((*process_countdown == 0 && this->create_shared_memory_)||
        (*process_countdown == 1 && !this->create_shared_memory_)){
      // If there is no process left using the shared memory or if the process
      // that created the shared memory is the only one use it then
      // delete the shared memory
      LOG(INFO) <<
          "This is the last process using shared memory so removing memory: " <<
          this->shared_memory_name_;
      bi::shared_memory_object::remove(
          this->shared_memory_name_.c_str());
    } else if (!this->create_shared_memory_){  // if this is not the process
      // that created the shared memory.
      // The countdown will NOT be reduced if the shared memory was created
      // by this process
      *process_countdown -= 1;
      LOG(INFO) << "This is not the last process using shared memory" <<
          "(remaining: " << *process_countdown <<") so keeping memory open: "<<
          this->shared_memory_name_;
    }
    delete this->shared_memory_;
    this->shared_memory_ = nullptr;
    this->shared_memory_name_ = "";
  }
}

GbzReader::GbzReader(const std::string& gbz_path,
                     const std::string& sample_name,
                     int context,
                     const std::string& chrom_prefix,
                     const std::string& shared_memory_name,
                     bool create_shared_memory,
                     bool use_loaded_shared_memory,
                     int shared_memory_size_gb,
                     int num_processes)
    : sample_name_(sample_name),
      context_(context),
      chrom_prefix_(chrom_prefix),
      shared_memory_name_(shared_memory_name),
      shared_memory_size_gb_(shared_memory_size_gb),
      create_shared_memory_(create_shared_memory),
      use_loaded_shared_memory_(use_loaded_shared_memory),
      num_processes_(num_processes){
  double start = gbwt::readTimer();

  // Open GBZ file in read mode.
  std::ifstream in(gbz_path);

  this->shared_memory_ = nullptr;
  bool shared_memory_is_used = (this->create_shared_memory_ ||
                                this->use_loaded_shared_memory_);

  // check if we are using shared memory or not
  if (shared_memory_is_used){
    // Create/Open a shared memory segment.
    this->create_or_open_shared_memory();
    // Create/Load a GBZ object in shared memory.
    this->gbz_shared_mem_ =
        std::make_unique<gbwtgraph::GBZ<gbwt::SharedMemCharAllocatorType>>(
            this->shared_memory_
        );
    // Load the GBZ file into the GBZ object.
    this->gbz_shared_mem_->simple_sds_load(in);
  } else {
    // Create a GBZ object in process memory.
    this->gbz_ = std::make_unique<gbwtgraph::GBZ<std::allocator<char>>>();
    // Load the GBZ file into the GBZ object.
    this->gbz_->simple_sds_load(in);
  }

  // Create a PathIndex object.
  this->path_index_ = shared_memory_is_used ?
      std::make_unique<gbwtgraph::PathIndex>(*this->gbz_shared_mem_) :
      std::make_unique<gbwtgraph::PathIndex>(*this->gbz_);

  if (shared_memory_ == nullptr && shared_memory_is_used){
    LOG(FATAL) << "shared memory could NOT be opened/created!\n";
  }

  if (this->create_shared_memory_ && this->shared_memory_ != nullptr){
    bool* stringarray_forward_only_loaded =
        shared_memory_->find_or_construct<bool>
        ("StringArray_forward_only_loaded")();
    bool* stringarray_final_loaded =
        shared_memory_->find_or_construct<bool>
        ("StringArray_final_loaded")();
    // We should set these two variables to true because it might not be
    // the first time that find_or_construct is called for these two bool
    // variables.
    // When we run load_gbz_into_shared_memory first it checks the existence of
    // these variables in check_existence_in_shared_memory() in
    // gbwt/src/support.cpp and creates these two variables there.
    *stringarray_forward_only_loaded = true;
    *stringarray_final_loaded = true;
    LOG(INFO) << "shared memory is created and strings are loaded\n";
  }

  double end = gbwt::readTimer();
  LOG(INFO) << "Loading GBZ file took " << end - start << " seconds";
}

nucleus::StatusOr<std::vector<nucleus::genomics::v1::Read>> GbzReader::Query(
    const Range& range) {

    const std::string& contig_name = range.reference_name();
    const size_t start_pos = range.start();
    const size_t end_pos = range.end();

    if (start_pos >= cache_start_pos_ + 300 &&
        end_pos <= std::max(cache_end_pos_ - 300, 0)) {
      return reads_cache_;
    }
    bool shared_memory_is_used = (this->use_loaded_shared_memory_
                                || this->create_shared_memory_);
    const gbwt::Metadata& metadata =
        shared_memory_is_used ? gbz_shared_mem_->index.metadata :
                                gbz_->index.metadata;

    // remove the prefix from the contig name
    std::string contig_name_without_prefix =
        contig_name.substr(chrom_prefix_.length());

    gbwt::size_type contig_id = metadata.contig(contig_name_without_prefix);
    // If the contig is not found, return an empty vector of reads
    if (contig_id == metadata.contig_names.size()) {
      std::vector<nucleus::genomics::v1::Read> reads_empty;
      return reads_empty;
    }
    // Get the path ids for the sample and contig
    std::vector<gbwt::size_type> path_ids = metadata.findPaths(
        metadata.sample(sample_name_),
        contig_id);

    if (path_ids.empty()) {
      return nucleus::NotFound(absl::StrCat(
          "Pangenome path ids not found for pangenome sample name: ",
          sample_name_));
    }
    handlegraph::path_handle_t path =
        shared_memory_is_used ?
        gbz_shared_mem_->graph.path_to_handle(path_ids.front()):
        gbz_->graph.path_to_handle(path_ids.front());

    gbwtgraph::SubgraphQuery query = gbwtgraph::SubgraphQuery::path_interval(
        path, start_pos, end_pos, context_,
        gbwtgraph::SubgraphQuery::HaplotypeOutput::all_haplotypes);


    std::vector<nucleus::genomics::v1::Read> reads;
    if (shared_memory_is_used){
      gbwtgraph::Subgraph subgraph(*gbz_shared_mem_, path_index_.get(), query);
      reads = GetReadsFromSubgraph(subgraph, *gbz_shared_mem_, chrom_prefix_);
    }else{
      gbwtgraph::Subgraph subgraph(*gbz_, path_index_.get(), query);
      reads = GetReadsFromSubgraph(subgraph, *gbz_, chrom_prefix_);
    }

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

template <typename CharAllocatorType>
std::vector<nucleus::genomics::v1::Read> GbzReader::GetReadsFromSubgraph(
    const gbwtgraph::Subgraph& subgraph,
    const gbwtgraph::GBZ<CharAllocatorType>& gbz,
    const std::string& chrom_prefix) {
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
        MakeRead(std::string(chrom_prefix) + contig_name,
                 start,
                 bases,
                 cigar_elements,
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

template <typename CharAllocatorType>
std::string GbzReader::GetBases(const gbwt::vector_type& path,
                                const gbwtgraph::GBZ<CharAllocatorType>& gbz) {
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
