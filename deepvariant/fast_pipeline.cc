/*
 * Copyright 2024 Google LLC.
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

#include "deepvariant/fast_pipeline.h"

#include <cstdlib>
#include <filesystem> // NOLINT
#include <fstream>
#include <ios>
#include <memory>
#include <string>


#include "deepvariant/fast_pipeline_utils.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "boost/asio.hpp" // NOLINT
#include "boost/interprocess/allocators/allocator.hpp" // NOLINT
#include "boost/interprocess/detail/os_file_functions.hpp" // NOLINT
#include "boost/interprocess/managed_shared_memory.hpp" // NOLINT
#include "boost/interprocess/shared_memory_object.hpp" // NOLINT
#include "boost/interprocess/sync/interprocess_semaphore.hpp" // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp" // NOLINT
#include "boost/interprocess/sync/scoped_lock.hpp" // NOLINT
#include "boost/process.hpp" // NOLINT
#include "boost/process/search_path.hpp" // NOLINT

ABSL_FLAG(std::string, make_example_flags, "",
          "file containing make_examples flags");
ABSL_FLAG(std::string, call_variants_flags, "",
          "file containing call_variants flags");
ABSL_FLAG(std::string, shm_prefix, "", "prefix for shared memory objects");
ABSL_FLAG(int, num_shards, 0, "number of make_examples shards");
ABSL_FLAG(int, buffer_size, 10485760,
          "Shared memory buffer size for each shard, default is 10MB");

namespace learning {
namespace genomics {
namespace deepvariant {

namespace bp = boost::process;
namespace bi = boost::interprocess;

FastPipeline::FastPipeline(int num_shards, int buffer_size,
                           absl::string_view shm_prefix,
                           absl::string_view path_to_make_examples_flags,
                           absl::string_view path_to_call_variants_flags,
                           absl::string_view dv_bin_path)
    : num_shards_(num_shards),
      buffer_size_(buffer_size),
      shm_prefix_(shm_prefix),
      dv_bin_path_(dv_bin_path) {
  shm_.resize(num_shards_);
  make_examples_processes_.resize(num_shards_);
  call_variants_processes_.resize(1);
  buffer_empty_.resize(num_shards_);
  items_available_.resize(num_shards_);
  make_examples_shard_finished_.resize(num_shards_);
  std::ifstream make_examples_flags_file;
  std::ifstream call_variants_flags_file;
  make_examples_flags_file.open(path_to_make_examples_flags.data(),
                                std::ios::in);
  call_variants_flags_file.open(path_to_call_variants_flags.data(),
                                std::ios::in);
  if (make_examples_flags_file.is_open()) {
    std::string line;
    while (std::getline(make_examples_flags_file, line)) {
      make_examples_flags_.push_back(line);
    }
  }
  if (call_variants_flags_file.is_open()) {
    std::string line;
    while (std::getline(call_variants_flags_file, line)) {
      call_variants_flags_.push_back(line);
    }
  }
}

// Set shared memory buffer and mutexes for all shards. These objects will be
// accessed by make_examples and call_variants.
void FastPipeline::SetGlobalObjects() {
  for (int shard = 0; shard < num_shards_; ++shard) {
    // Create shared memory buffers.
    std::string shard_shm_name = GetShmBufferName(shm_prefix_, shard);
    shm_[shard] =
        std::make_unique<bi::shared_memory_object>(bi::shared_memory_object(
            bi::open_or_create, shard_shm_name.data(), bi::read_write));
    shm_[shard]->truncate(buffer_size_);
    // Create mutex signalling that buffer is empty.
    LOG(INFO) << "Creating buffer_empty mutex";
    buffer_empty_[shard] = std::make_unique<bi::named_mutex>(
        bi::open_or_create,
        GetBufferEmptyMutexName(shm_prefix_, shard).data());
    // Create mutex signalling that items are available in the buffer.
    LOG(INFO) << "Creating items_available mutex";
    items_available_[shard] = std::make_unique<bi::named_mutex>(
        bi::open_or_create,
        GetItemsAvailableMutexName(shm_prefix_, shard).data());
    // Create mutex signalling that shard is finished.
    LOG(INFO) << "Creating shard_finished mutex";
    make_examples_shard_finished_[shard] = std::make_unique<bi::named_mutex>(
        bi::open_or_create,
        GetShardFinishedMutexName(shm_prefix_, shard).data());
  }
}

void FastPipeline::ClearGlobalObjects() {
  for (int shard = 0; shard < num_shards_; ++shard) {
    shm_[shard].release();
    buffer_empty_[shard].release();
    items_available_[shard].release();
    make_examples_shard_finished_[shard].release();
    shm_[shard]->remove(GetShmBufferName(shm_prefix_, shard).data());
    buffer_empty_[shard]->remove(
        GetBufferEmptyMutexName(shm_prefix_, shard).data());
    items_available_[shard]->remove(
        GetItemsAvailableMutexName(shm_prefix_, shard).data());
    make_examples_shard_finished_[shard]->remove(
        GetShardFinishedMutexName(shm_prefix_, shard).data());
  }
}

// Spawn make_examples processes.
void FastPipeline::SpawnMakeExamples() {
  std::filesystem::path make_examples_bin = dv_bin_path_;
  make_examples_bin /= "make_examples";
  LOG(INFO) << "make_examples_bin: " << make_examples_bin;
  make_examples_flags_.push_back(absl::StrCat("--shm_prefix=", shm_prefix_));
  make_examples_flags_.push_back(
      absl::StrCat("--shm_buffer_size=", buffer_size_));
  make_examples_flags_.push_back("--stream_examples");
  for (int shard = 0; shard < num_shards_; ++shard) {
    LOG(INFO) << "Spawning make_examples process for shard " << shard;
    make_examples_processes_[shard] =
        std::make_unique<bp::child>(bp::child(
            make_examples_bin.c_str(),
            make_examples_flags_,
            absl::StrCat("--task=", shard),
            bp::std_out > bp::null,
            bp::std_err > stderr));
    LOG(INFO) << "make_examples process " << shard << " process stared";
  }
}

// Spawn call_variants processes.
void FastPipeline::SpawnCallVariants() {
  std::filesystem::path call_variants_bin = dv_bin_path_;
  call_variants_bin /= "call_variants";
  LOG(INFO) << "call_variants_bin: " << call_variants_bin;
  call_variants_flags_.push_back(absl::StrCat("--shm_prefix=", shm_prefix_));
  call_variants_flags_.push_back("--stream_examples");
  call_variants_flags_.push_back(
      absl::StrCat("--num_input_shards=", num_shards_));
  LOG(INFO) << "Spawning call_variants process";
  call_variants_processes_[0] =
      std::make_unique<bp::child>(
          bp::child(call_variants_bin.c_str(),
                    call_variants_flags_,
                    bp::std_out > bp::null,
                    bp::std_err > stderr));
  LOG(INFO) << "call_variants process stared";
}

void FastPipeline::WaitForProcesses() {
  // TODO Maybe there is a better way to wait for all processes. Groups?
  for (int shard = 0; shard < num_shards_; ++shard) {
    make_examples_processes_[shard]->wait();
  }
  call_variants_processes_[0]->wait();
}

void RunFastPipeline(absl::string_view dv_bin_path) {
  FastPipeline fast_pipeline(
      absl::GetFlag(FLAGS_num_shards), absl::GetFlag(FLAGS_buffer_size),
      absl::GetFlag(FLAGS_shm_prefix), absl::GetFlag(FLAGS_make_example_flags),
      absl::GetFlag(FLAGS_call_variants_flags), dv_bin_path);
  fast_pipeline.SetGlobalObjects();
  fast_pipeline.SpawnMakeExamples();
  fast_pipeline.SpawnCallVariants();
  fast_pipeline.WaitForProcesses();
  fast_pipeline.ClearGlobalObjects();
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char** argv) {

  absl::ParseCommandLine(argc, argv);

  const char* dv_bin_path = std::getenv("DV_BIN_PATH");
  if (dv_bin_path == nullptr) {
    LOG(FATAL) << "DV_BIN_PATH environment variable is not set";
  }

  // TODO Add a check for required flags.
  // 1. shm_prefix has no spaces and match the pattern [a-zA-Z0-9_-]
  // 2. num_shards > 0
  // 3. buffer_size > 1M and < 1G (Check that buffer_size * num_shards can fit
  // into memory)
  // 4. make_example_flags file exists
  // 5. call_variants_flags file exists
  // 6. No SHM files with the same prefix exist.

  learning::genomics::deepvariant::RunFastPipeline(dv_bin_path);
  return EXIT_SUCCESS;
}
