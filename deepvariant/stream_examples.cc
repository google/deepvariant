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

#include "deepvariant/stream_examples.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "deepvariant/pileup_image_native.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Variant = nucleus::genomics::v1::Variant;

StreamExamples::StreamExamples(
    const MakeExamplesOptions& options,
    const AltAlignedPileup& alt_aligned_pileup)
    : options_(options), alt_aligned_pileup_(alt_aligned_pileup) {
  // Calculate the size of the pileup image.
  pileup_image_size_ = options_.pic_options().width() *
                       options_.pic_options().height() *
                       options_.pic_options().channels().size();
  if (alt_aligned_pileup_ == AltAlignedPileup::kRows) {
    pileup_image_size_ *= 3;
  }
  absl::string_view shm_prefix = options.shm_prefix();
  // Name of the shared memory buffer.
  shm_name_ = absl::StrCat(shm_prefix, "_shm_", options_.task_id());
  shm_buffer_size_ = options_.shm_buffer_size();

  // Open an existing shared memory buffer.
  shm_ = boost::interprocess::shared_memory_object(
      boost::interprocess::open_only, shm_name_.data(),
      boost::interprocess::read_write);
  // Map the shared memory buffer to the address space of the process.
  shm_region_ =
      boost::interprocess::mapped_region(shm_, boost::interprocess::read_write);
  shm_buffer_ = static_cast<unsigned char*>(shm_region_.get_address());
  shm_buffer_pos_ = 0;

  // Init shared mutexes.
  LOG(INFO) << "Creating buffer_empty mutex";
  buffer_empty_ = std::make_unique<boost::interprocess::named_mutex>(
      boost::interprocess::open_only,
      absl::StrCat(shm_prefix, "_buffer_empty_", options_.task_id()).data());
  LOG(INFO) << "Creating items_available mutex";
  items_available_ = std::make_unique<boost::interprocess::named_mutex>(
      boost::interprocess::open_only,
      absl::StrCat(shm_prefix, "_items_available_", options_.task_id()).data());
  items_available_->lock();
  LOG(INFO) << "Creating shard_finished mutex";
  shard_finished_ = std::make_unique<boost::interprocess::named_mutex>(
      boost::interprocess::open_only,
      absl::StrCat(shm_prefix, "_shard_finished_", options_.task_id()).data());
  shard_finished_->lock();
}

// Writes len as a sequence of bytes to the shm_buffer_.
void StreamExamples::WriteLenToShm(int len) {
  std::copy(
      static_cast<const char*>(static_cast<const void*>(&len)),
      static_cast<const char*>(static_cast<const void*>(&len)) + sizeof len,
      &shm_buffer_[0] + shm_buffer_pos_);
  shm_buffer_pos_ += sizeof len;
}

// Writes the length of the data first, then writes `len` number of bytes to the
// shm_buffer_.
void StreamExamples::WriteBytesToShm(const void* bytes, int len) {
  WriteLenToShm(len);
  std::copy(static_cast<const char*>(bytes),
            static_cast<const char*>(bytes) + len,
            &shm_buffer_[0] + shm_buffer_pos_);
  shm_buffer_pos_ += len;
}

// Puts the example into the shared memory buffer. Once the buffer is full,
// waits until call_variants is done processing the buffer and then continues.
// The recipient is notified once the buffer is full or all examples are
// processed.
void StreamExamples::StreamExample(
    const std::vector<std::unique_ptr<ImageRow>>& ref_rows,
    const std::vector<std::vector<std::unique_ptr<ImageRow>>>& alt_images,
    const AltAlignedPileup& alt_aligned_pileup,
    absl::string_view alt_indices_encoded, absl::string_view variant_encoded) {
  // Calculate the space needed.
  int required_space = variant_encoded.size() + alt_indices_encoded.size() +
                       pileup_image_size_ + (sizeof(int) * 4);

  while (true) {
    // Check is there is enough space left in the shm_buffer_ to add the
    // example.
    if (required_space < shm_buffer_size_ - shm_buffer_pos_) {
      // Write alt_indices serialized.
      WriteBytesToShm(alt_indices_encoded.data(), alt_indices_encoded.size());
      // Write variant serialized.
      WriteBytesToShm(variant_encoded.data(), variant_encoded.size());
      // Write pileup image size.
      WriteLenToShm(pileup_image_size_);
      // Write pileup image. FillPileupArray populates the shm_buffer_ directly.
      FillPileupArray(ref_rows, alt_images, alt_aligned_pileup, &shm_buffer_,
                      shm_buffer_size_, shm_buffer_pos_);
      shm_buffer_pos_ += pileup_image_size_;
      break;
    } else {
      // Let call_variants know that items are ready and wait for the buffer to
      // become available.
      if (shm_buffer_size_ - shm_buffer_pos_ > 0) {
        // Write 0 to indicate the end of the batch.
        WriteLenToShm(0);
      }
      items_available_->unlock();
      // Wait for buffer to become empty.
      buffer_empty_->lock();
      shm_buffer_pos_ = 0;
    }
  }
}

void StreamExamples::StartStreaming() {
  buffer_empty_->lock();
  shm_buffer_pos_ = 0;
}

void StreamExamples::EndStreaming(bool data_written) {
  if (data_written) {
    WriteLenToShm(0);
    items_available_->unlock();
  } else {
    buffer_empty_->unlock();
  }
}

void StreamExamples::SignalShardFinished() {
    buffer_empty_->lock();
    shard_finished_->unlock();
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
