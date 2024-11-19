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

#ifndef THIRD_PARTY_NUCLEUS_IO_GBZ_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_GBZ_READER_H_

#include <memory>
#include <vector>
#include <string>

#include "include/gbwt/utils.h"
#include "include/gbwtgraph/gbz.h"
#include "include/gbwtgraph/subgraph.h"
#include "third_party/nucleus/core/statusor.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

#include "boost/interprocess/managed_shared_memory.hpp"  // NOLINT
#include "boost/interprocess/shared_memory_object.hpp"  // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp"  // NOLINT


namespace nucleus {

using SamIterable = Iterable<nucleus::genomics::v1::Read>;

// A GBZ reader.
//
// GBZ files store information about pangenome population data.
//
// This class provides methods to iterate through a GBZ file and
// query() for only read overlapping a specific region on the genome.
//
// The objects returned by iterate() or query() are nucleus.genomics.v1.Read
// objects parsed from the GBZ records in the file.
class GbzReader : public Reader {
 public:
  // Creates a new GbzReader reading from the GBZ file gbz_path.
  //
  // gbz_path must point to an existing GBZ formatted file.
  //
  // Returns a StatusOr that is OK if the GbzReader could be successfully
  // created or an error code indicating the error that occurred.
  GbzReader(const std::string& gbz_path,
            const std::string& sample_name,
            int context,
            const std::string& chrom_prefix = "",
            const std::string& shared_memory_name = "GBZ_SHARED_MEMORY",
            bool create_shared_memory = false,
            bool use_loaded_shared_memory = false,
            int shared_memory_size_gb = 12,
            int num_processes = 0);

  ~GbzReader(){this->close_shared_memory();}

  nucleus::StatusOr<std::vector<nucleus::genomics::v1::Read>> Query(
      const nucleus::genomics::v1::Range& range);

  nucleus::StatusOr<std::shared_ptr<SamIterable>> Iterate() const;


 private:
  // The filename of the GBZ file.
  std::string gbz_filename_;
  // The sample name of the sample for query.
  std::string sample_name_;

  // Below are the pointers to the GBZ objects in shared memory and in process
  // memory. Only one of them will be used at a time.
  // when shared memory is used, gbz_shared_mem_ will be used.
  // when shared memory is not used, gbz_ will be used.
  // (TODO mobinasri: maybe use template class later to keep only one pointer)
  // Pinter to the GBZ object in shared memory.
  std::unique_ptr<gbwtgraph::GBZ<gbwt::SharedMemCharAllocatorType>>
      gbz_shared_mem_;
  // Pinter to the GBZ object in process memory.
  std::unique_ptr<gbwtgraph::GBZ<std::allocator<char>>>
      gbz_;

  // The PathIndex  object.
  std::unique_ptr<gbwtgraph::PathIndex> path_index_;
  // context size
  int context_;
  // chrom prefix
  std::string chrom_prefix_;

  // shared memory
  boost::interprocess::managed_shared_memory* shared_memory_;
  // shared memory name
  std::string shared_memory_name_;
  // shared memory size
  int shared_memory_size_gb_;
  // If True, create shared memory
  bool create_shared_memory_;
  // If True, use shared memory that is created by other process
  bool use_loaded_shared_memory_;
  // number of processes that will use shared memory excluding the one that
  // created it.
  // it is used to make sure that the shared memory is not deleted before all
  // processes are done using it.
  int num_processes_;

  void create_or_open_shared_memory();
  void close_shared_memory();

  std::vector<nucleus::genomics::v1::Read> reads_cache_;
  int cache_start_pos_ = 0;
  int cache_end_pos_ = 0;

  void updateCache(const std::vector<nucleus::genomics::v1::Read>& reads);

  template <typename CharAllocatorType>
  static std::string GetBases(const gbwt::vector_type& path,
                              const gbwtgraph::GBZ<CharAllocatorType>& gbz);
  static std::string GetReverseComplement(const std::string& sequence);
  static nucleus::genomics::v1::Read MakeRead(
      const std::string& chr, const int start, const std::string& bases,
      const std::vector<std::string>& cigar_elements,
      const std::string& read_name);
  static genomics::v1::CigarUnit_Operation ParseCigarOpStr(const char op);
  static std::vector<std::string> GetCigarElements(const std::string& input);

  template <typename CharAllocatorType>
  static std::vector<nucleus::genomics::v1::Read> GetReadsFromSubgraph(
      const gbwtgraph::Subgraph& subgraph,
      const gbwtgraph::GBZ<CharAllocatorType>& gbz,
      const std::string& chrom_prefix);
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_GBZ_READER_H_
