/*
 * Copyright 2022 Google LLC.
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

#ifndef THIRD_PARTY_NUCLEUS_IO_VARIANT_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_VARIANT_READER_H_

#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "third_party/nucleus/io/tfrecord_reader.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace nucleus {

constexpr std::string_view kAutoDetectCompression = "AUTO";

using nucleus::genomics::v1::Variant;

// Holds a pointer to the Variant proto, and the index the contig it belongs to.
struct IndexedVariant {
  std::unique_ptr<Variant> variant;
  uint32_t contig_map_index;

  bool operator>(const IndexedVariant& other) const;
};

IndexedVariant EmptyIndexedVariant();

// Reads Variant proto records from a single TFRecord file.
//
// The index of the contig each variant belongs to is returned together with it,
// in order to simplify contig-index based ordering later on.
//
// Note: This is intended for a specific use case within DeepVariant, where both
// the Variant and the index of the contig are always used together.
// If you need a more general use case, consider using TFRecordReader directly.
class VariantReader {
 public:
  // Internal constructor, `Open` should generally be used instead.
  VariantReader(std::unique_ptr<TFRecordReader> internal_reader,
                absl::flat_hash_map<std::string, uint32_t>& contig_index_map);

  // Creates a reader for the given file.
  // `compression_type` can be either "" (for no compression), "GZIP", or "AUTO"
  // (for auto detection by filename suffix).
  // `contig_index_map` should be a mapping between Variant reference names and
  // their index within the sorted contigs.
  static std::unique_ptr<VariantReader> Open(
      const std::string& filename, std::string_view compression_type,
      absl::flat_hash_map<std::string, uint32_t>& contig_index_map);

  IndexedVariant GetAndReadNext();

  // Reads the next record if available.
  bool GetNext();

  // Returns the current Variant and contig index.
  // Only valid after GetNext() has returned true.
  IndexedVariant ReadRecord();

 private:
  std::unique_ptr<TFRecordReader> internal_reader_;
  absl::flat_hash_map<std::string, uint32_t> contig_index_map_;
};

struct VariantFromShard {
  // This allows accessing the inner unique_ptr to the variant even when stored
  // in std::priority_queue, which sadly only supports const ref .top() .
  // But note that if that unique_ptr is moved, the element must be poped from
  // the priority_queue right after that (or the heap will try to compare
  // invalid data).
  mutable IndexedVariant variant;
  uint32_t reader_shard_index;
};

// Ranking function for priority_queue. Using a > b allows it to act as a
// min_heap and not like a max_heap as it would by default.
struct CompareVariantFromShard {
  bool operator()(const VariantFromShard& a, const VariantFromShard& b) const {
    return a.variant > b.variant;
  }
};

// Reads Variant proto records from sharded TFRecord file paths in sorted order.
//
// The input TFRecord file must have each shard already in sorted order (but
// elements can be interleaved across shards). Under those constraints, the
// elements will be returned in a global sorted order.
//
// Note: This is intended for a specific use case within DeepVariant, where both
// the Variant and the index of the contig are always used together and are used
// for sorting the Variants.
// For a more general use case, consider using multiple TFRecordReaders
// directly.
class ShardedVariantReader {
 public:
  // Internal constructor, `Open` should generally be used instead.
  ShardedVariantReader(
      std::vector<std::unique_ptr<VariantReader>> shard_readers);

  // Creates a reader for the given file paths.
  // `compression_type` can be either "" (for no compression), "GZIP", or "AUTO"
  // (for auto detection by filename suffix). `contig_index_map` should be a
  // mapping between reference names and their index within the sorted contigs.
  static std::unique_ptr<ShardedVariantReader> Open(
      const std::vector<std::string>& shard_paths,
      absl::flat_hash_map<std::string, uint32_t>& contig_index_map);

  IndexedVariant GetAndReadNext();

 private:
  bool GetNext();
  void ReadNextFromShard(uint32_t shard_idx);

  IndexedVariant next_elem_;
  // Min_heap which yields the next *globally* 'smallest' Variant each time.
  std::priority_queue<VariantFromShard, std::vector<VariantFromShard>,
                      CompareVariantFromShard>
      next_elems_;
  std::vector<std::unique_ptr<VariantReader>> shard_readers_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_VARIANT_READER_H_
