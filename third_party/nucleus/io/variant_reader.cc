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

#include "third_party/nucleus/io/variant_reader.h"

#include <cstdint>
#include <memory>

#include "absl/log/check.h"
#include "absl/strings/match.h"
#include "third_party/nucleus/io/tfrecord_reader.h"

namespace nucleus {

bool IndexedVariant::operator>(const IndexedVariant& other) const {
  if (contig_map_index > other.contig_map_index) {
    return true;
  }
  if (contig_map_index < other.contig_map_index) {
    return false;
  }
  if (variant->start() > other.variant->start()) {
    return true;
  }
  if (variant->start() < other.variant->start()) {
    return false;
  }
  return variant->end() > other.variant->end();
}

IndexedVariant EmptyIndexedVariant() {
  return {.variant = nullptr,
          // An empty empty value is comes 'later' (sorting wise) than any other
          // value, so its index should be infinity.
          .contig_map_index = std::numeric_limits<uint32_t>::max()};
}

VariantReader::VariantReader(
    std::unique_ptr<TFRecordReader> internal_reader,
    absl::flat_hash_map<std::string, uint32_t>& contig_index_map)
    : internal_reader_(std::move(internal_reader)),
      contig_index_map_(contig_index_map) {}

std::unique_ptr<VariantReader> VariantReader::Open(
    const std::string& filename, std::string_view compression_type,
    absl::flat_hash_map<std::string, uint32_t>& contig_index_map) {
  std::string compression(compression_type);
  if (compression_type == kAutoDetectCompression) {
    compression = "";
    if (absl::EndsWith(filename, ".gz")) {
      compression = "GZIP";
    }
  }

  return std::make_unique<VariantReader>(
      TFRecordReader::New(filename, compression), contig_index_map);
}

bool VariantReader::GetNext() { return internal_reader_->GetNext(); }

// Return the current record contents.  Only valid after GetNext()
// has returned true.
IndexedVariant VariantReader::ReadRecord() {
  tensorflow::tstring data = internal_reader_->record();
  std::unique_ptr<Variant> proto = std::make_unique<Variant>();
  CHECK(proto->ParseFromArray(data.data(), data.length()))
      << "Failed to parse proto";
  uint32_t contig_index = contig_index_map_[proto->reference_name()];
  return {.variant = std::move(proto), .contig_map_index = contig_index};
}

IndexedVariant VariantReader::GetAndReadNext() {
  if (GetNext()) {
    return ReadRecord();
  }
  return EmptyIndexedVariant();
}

ShardedVariantReader::ShardedVariantReader(
    std::vector<std::unique_ptr<VariantReader>> shard_readers)
    : shard_readers_(std::move(shard_readers)) {
  // Prime all readers and 'next_elems_'
  for (uint32_t i = 0; i < shard_readers_.size(); i++) {
    ReadNextFromShard(i);
  }
}

std::unique_ptr<ShardedVariantReader> ShardedVariantReader::Open(
    const std::vector<std::string>& shard_paths,
    absl::flat_hash_map<std::string, uint32_t>& contig_index_map) {
  std::vector<std::unique_ptr<VariantReader>> shard_readers;
  shard_readers.reserve(shard_paths.size());
  for (const auto& path : shard_paths) {
    shard_readers.emplace_back(
        VariantReader::Open(path, kAutoDetectCompression, contig_index_map));
  }

  return std::make_unique<ShardedVariantReader>(std::move(shard_readers));
}

// All shards are sorted internally, but not with respect to each other.
bool ShardedVariantReader::GetNext() {
  if (next_elems_.empty()) {
    return false;
  }

  const auto& next = next_elems_.top();
  int min_elem_idx = next.reader_shard_index;
  next_elem_ = std::move(next.variant);
  next_elems_.pop();
  ReadNextFromShard(min_elem_idx);

  return true;
}

IndexedVariant ShardedVariantReader::GetAndReadNext() {
  if (GetNext()) {
    return std::move(next_elem_);
  }
  return EmptyIndexedVariant();
}

void ShardedVariantReader::ReadNextFromShard(uint32_t shard_idx) {
  // Advance the reader which was used, and remove if fully consumed
  if (shard_readers_[shard_idx]->GetNext()) {
    IndexedVariant variant = shard_readers_[shard_idx]->ReadRecord();
    next_elems_.push(
        {.variant = std::move(variant), .reader_shard_index = shard_idx});
  }
}

}  // namespace nucleus
