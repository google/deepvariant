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

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#define NO_SANITIZE __attribute__((no_sanitize("all")))

#include "absl/container/flat_hash_set.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/synchronization/mutex.h"
#include "boost/interprocess/managed_shared_memory.hpp" // NOLINT
#include "boost/interprocess/shared_memory_object.hpp" // NOLINT
#include "boost/interprocess/sync/named_mutex.hpp" // NOLINT
#include "tensorflow/core/framework/op_kernel.h"
#include "tensorflow/core/framework/op_requires.h"
#include "tensorflow/core/framework/resource_mgr.h"  // copybara: quotes ""
#include "tensorflow/core/framework/resource_op_kernel.h"  // copybara: quotes ""
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/framework/tensor_shape.h"
#include "tensorflow/core/framework/types.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/mutex.h"
#include "tensorflow/core/platform/refcount.h"
#include "tensorflow/core/platform/status.h"
#include "tensorflow/core/platform/tstring.h"
#include "tensorflow/core/platform/types.h"
#include "tensorflow/tsl/platform/errors.h"
#include "tensorflow/tsl/platform/thread_annotations.h"
#include "deepvariant/fast_pipeline_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// This class implements a custom TensorFlow op that reads data from shared
// memory buffers and returns a batch of examples. Batches are variable length.
// Examples are read from shared memory buffers (one for each shard).
// Examples are written in the following format:
//   int: length of alt_indices_encoded
//   string: alt_indices_encoded
//   int: length of variant_encoded
//   string: variant_encoded
//   int: length of raw pileup_image
//   string: array containing pileup image.
// The operation is utilized by call_variants input dataset. Dataset
// parallelizes extraction using multiple threads. This way multiple threads
// make calls to operation `Next` in parallel. Each instance of Next tris to
// lock the mutex for the next available shard. If mutex cannot be acquired
// then the next shard is tried. If all shards are completed then the operation
// returns an empty batch.
// TODO Add unit tests for Init and Next operatios.
class StreamExamplesResource : public tensorflow::ResourceBase {
 public:
  explicit StreamExamplesResource(tensorflow::Env* env) : env_(env) {}

  virtual tensorflow::Status Init(const tensorflow::string& shm_prefix,
                                  int num_shards) {
    tensorflow::mutex_lock l(mu_);
    shm_buffer_.resize(num_shards);
    shm_.resize(num_shards);
    shm_region_.resize(num_shards);
    buffer_empty_.resize(num_shards);
    items_available_.resize(num_shards);
    make_examples_shard_finished_.resize(num_shards);
    num_shards_ = num_shards;
    for (int shard = 0; shard < num_shards; shard++) {
      shm_[shard] = std::make_unique<boost::interprocess::shared_memory_object>(
          boost::interprocess::shared_memory_object(
              boost::interprocess::open_only,
              GetShmBufferName(shm_prefix, shard).data(),
              boost::interprocess::read_write));
      shm_region_[shard] = std::make_unique<boost::interprocess::mapped_region>(
          *shm_[shard], boost::interprocess::read_write);
      shm_buffer_[shard] =
          static_cast<unsigned char*>(shm_region_[shard]->get_address());

      // Init mutexes
      buffer_empty_[shard] = std::make_unique<boost::interprocess::named_mutex>(
          boost::interprocess::open_only,
          GetBufferEmptyMutexName(shm_prefix, shard).data());
      items_available_[shard] =
          std::make_unique<boost::interprocess::named_mutex>(
              boost::interprocess::open_only,
              GetItemsAvailableMutexName(shm_prefix, shard).data());
      make_examples_shard_finished_[shard] =
          std::make_unique<boost::interprocess::named_mutex>(
              boost::interprocess::open_only,
              GetShardFinishedMutexName(shm_prefix, shard).data());
    }

    return tensorflow::Status();
  }

  NO_SANITIZE
  int ReadIntFromBuffer(const unsigned char* byte_array, int pos, int len) {
    return *reinterpret_cast<const int*>(&byte_array[pos]);
  }

  NO_SANITIZE
  std::string ReadStringFromBuffer(const unsigned char* byte_array, int pos,
                                   int len) {
    return std::string(&byte_array[pos], &byte_array[pos] + len);
  }

  NO_SANITIZE
  tensorflow::Status Next(
      const int64_t index,
      std::function<tensorflow::Status(
          const tensorflow::TensorShape& shape, tensorflow::Tensor** image,
          tensorflow::Tensor** variant, tensorflow::Tensor** alt_allele_idx)>
          allocate_func) {
    int shard = index % num_shards_;
    std::vector<tensorflow::string> variant_records;
    std::vector<tensorflow::string> alt_allele_idx_records;
    std::vector<tensorflow::string> image_records;
    bool shard_is_completed = false;
    int num_shards_completed = 0;
    {
      absl::MutexLock lock(&mu_wip_);
      num_shards_completed = completed_shards_.size();
    }

    // We try to read data from the `index` shard first. If data is not
    // available we immediately move to the next shard as long as we still have
    // unfinished shards.
    while (num_shards_completed < num_shards_) {
      {
        absl::MutexLock lock(&mu_wip_);
        num_shards_completed = completed_shards_.size();
      }
      // Try to acquire items_available mutex without waiting. If it is not
      // available move to the next shard.
      if (items_available_[shard]->try_lock()) {
        // Read alt indices.
        int shm_buffer_pos = 0;
        int len = ReadIntFromBuffer(shm_buffer_[shard], shm_buffer_pos,
                                    sizeof(int));
        shm_buffer_pos += sizeof(int);
        while (len > 0) {
          // Read the length of alt_indices_encoded.
          std::string alt_indices_encoded =
              ReadStringFromBuffer(shm_buffer_[shard], shm_buffer_pos, len);
          alt_allele_idx_records.push_back(alt_indices_encoded);
          shm_buffer_pos += len;

          // Read variant.
          len = ReadIntFromBuffer(shm_buffer_[shard], shm_buffer_pos,
                                  sizeof(int));
          shm_buffer_pos += sizeof(int);
          std::string varint_encoded =
              ReadStringFromBuffer(shm_buffer_[shard], shm_buffer_pos, len);
          shm_buffer_pos += len;
          variant_records.push_back(varint_encoded);

          // Read image.
          len = ReadIntFromBuffer(shm_buffer_[shard], shm_buffer_pos,
                                  sizeof(int));
          shm_buffer_pos += sizeof(int);
          std::string image_bytes(&(shm_buffer_[shard])[shm_buffer_pos],
                                  &(shm_buffer_[shard])[shm_buffer_pos + len]);
          image_records.push_back(image_bytes);
          shm_buffer_pos += len;

          // Reading next record
          len = ReadIntFromBuffer(shm_buffer_[shard], shm_buffer_pos,
                                  sizeof(int));
          shm_buffer_pos += sizeof(int);
        }
        buffer_empty_[shard]->unlock();
        break;
      } else {
        // Check if shard is completed. In that we need to
        if (make_examples_shard_finished_[shard]->try_lock()) {
          LOG(INFO) << "Shard " << shard << " is completed";
          std::pair<absl::flat_hash_set<int>::iterator, bool> res;
          {
            absl::MutexLock lock(&mu_wip_);
            res = completed_shards_.insert(shard);
            num_shards_completed = completed_shards_.size();
          }
          shard_is_completed = true;
          break;
        }
        shard = (shard + 1) % num_shards_;
      }
    }
    // Add records as strings to the output

    // If all shards are processed we need to signal the end of the stream by
    // sending an empty tensor, by way of using empty records.
    // If we fell off from the while loop that means the shard is completed.
    tensorflow::TensorShape shape({static_cast<int64_t>(image_records.size())});
    tensorflow::Tensor* image_tensor;
    tensorflow::Tensor* variant_tensor;
    tensorflow::Tensor* alt_allele_idx_tensor;
    TF_RETURN_IF_ERROR(allocate_func(shape, &image_tensor, &variant_tensor,
                                     &alt_allele_idx_tensor));
    for (size_t i = 0; i < image_records.size(); i++) {
      variant_tensor->flat<tensorflow::tstring>()(i) = variant_records[i];
      alt_allele_idx_tensor->flat<tensorflow::tstring>()(i) =
          alt_allele_idx_records[i];
      image_tensor->flat<tensorflow::tstring>()(i) = image_records[i];
    }
    return tensorflow::Status();
  }

  tensorflow::string DebugString() const override {
    return "StreamFromShmFilesResource";
  }

 protected:
  mutable tensorflow::mutex mu_;
  tensorflow::Env* env_ TF_GUARDED_BY(mu_);
  absl::Mutex mu_wip_;
  std::vector<std::unique_ptr<boost::interprocess::shared_memory_object>> shm_;
  std::vector<std::unique_ptr<boost::interprocess::mapped_region>> shm_region_;
  std::vector<unsigned char*> shm_buffer_;
  std::vector<std::unique_ptr<boost::interprocess::named_mutex>> buffer_empty_;
  std::vector<std::unique_ptr<boost::interprocess::named_mutex>>
      items_available_;
  std::vector<std::unique_ptr<boost::interprocess::named_mutex>>
      make_examples_shard_finished_;
  int num_shards_;
  absl::flat_hash_set<int> completed_shards_;
};

class StreamExamplesInitOp
    : public tensorflow::ResourceOpKernel<StreamExamplesResource> {
 public:
  explicit StreamExamplesInitOp(tensorflow::OpKernelConstruction* context)
      : ResourceOpKernel<StreamExamplesResource>(context) {
    env_ = context->env();
  }

 private:
  void Compute(tensorflow::OpKernelContext* context) override {
    ResourceOpKernel<StreamExamplesResource>::Compute(context);

    const tensorflow::Tensor* shm_dir_tensor;
    OP_REQUIRES_OK(context, context->input("shm_dir", &shm_dir_tensor));
    const std::string& shm_dir =
        shm_dir_tensor->scalar<tensorflow::tstring>()();

    const tensorflow::Tensor* num_shards_tensor;
    OP_REQUIRES_OK(context, context->input("num_shards", &num_shards_tensor));
    int num_shards = num_shards_tensor->scalar<tensorflow::int32>()();

    OP_REQUIRES_OK(context, get_resource()->Init(shm_dir, num_shards));
  }
  tensorflow::Status CreateResource(StreamExamplesResource** resource)
      TF_EXCLUSIVE_LOCKS_REQUIRED(mu_) override {
    *resource = new StreamExamplesResource(env_);
    return tensorflow::Status();
  }

 private:
  // mu_ is derived from StreamExamplesResource.
  // mutable tensorflow::mutex mu_;
  tensorflow::Env* env_ TF_GUARDED_BY(mu_);
};

class StreamExamplesNextOp : public tensorflow::OpKernel {
 public:
  explicit StreamExamplesNextOp(tensorflow::OpKernelConstruction* context)
      : OpKernel(context) {
    env_ = context->env();
  }

  void Compute(tensorflow::OpKernelContext* context) override {
    StreamExamplesResource* resource;
    OP_REQUIRES_OK(context,
                   GetResourceFromContext(context, "input", &resource));
    tensorflow::core::ScopedUnref unref(resource);

    const tensorflow::Tensor* index_tensor;
    OP_REQUIRES_OK(context, context->input("index", &index_tensor));
    const int64_t index = index_tensor->scalar<int64_t>()();

    OP_REQUIRES_OK(
        context,
        resource->Next(
            index,
            [&](const tensorflow::TensorShape& shape,
                tensorflow::Tensor** image, tensorflow::Tensor** variant,
                tensorflow::Tensor** alt_allele_idx) -> tensorflow::Status {
              TF_RETURN_IF_ERROR(context->allocate_output(0, shape, image));
              TF_RETURN_IF_ERROR(context->allocate_output(1, shape, variant));
              TF_RETURN_IF_ERROR(
                  context->allocate_output(2, shape, alt_allele_idx));
              return tensorflow::Status();
            }));
  }

 private:
  mutable tensorflow::mutex mu_;
  tensorflow::Env* env_ TF_GUARDED_BY(mu_);
};

REGISTER_KERNEL_BUILDER(
    Name("StreamExamplesInit").Device(tensorflow::DEVICE_CPU),
    StreamExamplesInitOp);
REGISTER_KERNEL_BUILDER(
    Name("StreamExamplesNext").Device(tensorflow::DEVICE_CPU),
    StreamExamplesNextOp);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
