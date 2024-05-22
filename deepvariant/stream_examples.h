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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_STREAM_EXAMPLES_H_
#define LEARNING_GENOMICS_DEEPVARIANT_STREAM_EXAMPLES_H_

#include "deepvariant/pileup_image_native.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "boost/interprocess/managed_shared_memory.hpp" // NOLINT
#include "boost/interprocess/shared_memory_object.hpp" // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp" // NOLINT

namespace learning {
namespace genomics {
namespace deepvariant {

class StreamExamples {
 public:
  explicit StreamExamples(
      const MakeExamplesOptions& options,
      const AltAlignedPileup& alt_aligned_pileup);

  // Writes examples to shared memory buffer.
  void StreamExample(
      const std::vector<std::unique_ptr<ImageRow>>& ref_rows,
      const std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_images,
      const AltAlignedPileup& alt_aligned_pileup,
      absl::string_view alt_indices_encoded, absl::string_view variant_encoded);

  void StartStreaming();
  void EndStreaming(bool data_written);
  void SignalShardFinished();

  // Helper methods to write bytes to shared memory.
  void WriteBytesToShm(const void* bytes, int len);

  // Helper method to write length of the data to shared memory.
  void WriteLenToShm(int len);

 private:
  // Make examples config.
  const MakeExamplesOptions options_;

  // Alt aligned pileup option.
  AltAlignedPileup alt_aligned_pileup_;

  // Shared memory object.
  boost::interprocess::shared_memory_object shm_;

  // Shared memory object name (shared between all processes).
  std::string shm_name_;

  // Shared memory buffer size.
  int shm_buffer_size_;

  // Size of the pileup image array.
  int pileup_image_size_;

  // Shared memory buffer mapped to the current process memory.
  unsigned char* shm_buffer_;

  // Mapped region of the shared memory.
  boost::interprocess::mapped_region shm_region_;

  // Current position in the shared memory buffer.
  int shm_buffer_pos_;

  // Synchronization mutexes for streaming examples.
  std::unique_ptr<boost::interprocess::named_mutex> buffer_empty_;
  std::unique_ptr<boost::interprocess::named_mutex> items_available_;
  std::unique_ptr<boost::interprocess::named_mutex> shard_finished_;
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_STREAM_EXAMPLES_H_
