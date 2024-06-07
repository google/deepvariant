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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_FAST_PIPELINE_H_
#define LEARNING_GENOMICS_DEEPVARIANT_FAST_PIPELINE_H_

#include <memory>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "boost/interprocess/shared_memory_object.hpp" // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp" // NOLINT
#include "boost/process.hpp" // NOLINT

namespace learning {
namespace genomics {
namespace deepvariant {

class FastPipeline {
 public:
  FastPipeline(int num_shards, int buffer_size, absl::string_view shm_prefix,
               absl::string_view path_to_make_examples_flags,
               absl::string_view path_to_call_variants_flags,
               absl::string_view dv_bin_path);
  void SetGlobalObjects();
  void SpawnMakeExamples();
  void SpawnCallVariants();
  void WaitForProcesses();
  void ClearGlobalObjects();

 private:
  int num_shards_;
  int buffer_size_;
  absl::string_view shm_prefix_;
  absl::string_view dv_bin_path_;
  std::vector<std::unique_ptr<boost::interprocess::shared_memory_object>> shm_;

  std::vector<std::unique_ptr<boost::interprocess::named_mutex>> buffer_empty_;
  std::vector<std::unique_ptr<boost::interprocess::named_mutex>>
      items_available_;
  std::vector<std::unique_ptr<boost::interprocess::named_mutex>>
      make_examples_shard_finished_;

  std::vector<std::unique_ptr<boost::process::child>> make_examples_processes_;
  std::vector<std::unique_ptr<boost::process::child>> call_variants_processes_;

  std::vector<std::string> make_examples_flags_;
  std::vector<std::string> call_variants_flags_;
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_FAST_PIPELINE_H_
