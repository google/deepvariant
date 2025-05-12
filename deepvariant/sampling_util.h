/*
 * Copyright 2025 Google LLC.
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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_SAMPLING_UTIL_H_
#define LEARNING_GENOMICS_DEEPVARIANT_SAMPLING_UTIL_H_

#include <cstddef>
#include <functional>
#include <random>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"

namespace learning::genomics::deepvariant::sampling {

namespace internal {

//// Implementation of InPlaceReservoirSample, exposing the input randomness.
// `index_provider` is a function that takes an index into the population vector
// and returns the index of the element that should be swapped with the element
// at that index.
template <typename T>
int InPlaceReservoirSampleImpl(int sample_size,
                               std::function<size_t(size_t)> index_provider,
                               std::vector<T>& population) {
  // If the vector is already smaller than the sample size, return the size.
  if (population.size() < sample_size) {
    return population.size();
  }
  for (size_t index = sample_size; index < population.size(); ++index) {
    size_t swap_index = index_provider(index);
    if (swap_index < sample_size) {
      std::swap(population[swap_index], population[index]);
    }
  }
  return sample_size;
}

}  // namespace internal

//// Attempts to sample sample_size elements from the population without
/// replacement. If sample_size is larger than the population, the
/// entire population is returned. Sampling is done in place, as a prefix of the
/// population vector, and each possible sample is chosen uniformly at random.
template <typename T>
int InPlaceReservoirSample(int sample_size, std::mt19937_64& gen,
                           std::vector<T>& population) {
  auto uniform_index_dist = [&gen](size_t max) {
    auto dist = std::uniform_int_distribution<size_t>(0, max);
    return dist(gen);
  };
  return internal::InPlaceReservoirSampleImpl(sample_size, uniform_index_dist,
                                              population);
}

//// Samples from a family of vectors, where each vector represents a partition.
/// The function samples sample_size elements from the union of the vectors,
/// attempting to ensure that each partition has at least min_per_partition
/// elements in the sample.
///
/// NOTE: This algorithm doesn't produce a uniform distribution; it tends to
/// produce samples that are more `balanced`. For example, consider the case
/// where we have two partitions, each with 3 elements, and we want to sample 4
/// elements, with a minimum of 1 element per partition. There are 3^2 = 9
/// balanced possible choices, and 6 unbalanced ones, so if we were picking
/// uniformly, we would expect 3/5 of the samples to be balanced. In contrast,
/// this function will pick a single element from each partition first, and then
/// choose 2 elements from the union of the remaining 4 elements. This choice
/// will be balanced 2/3 of the time.
template <typename T>
absl::StatusOr<std::vector<T>> SampleWithPartitionMins(
    std::vector<std::vector<T>> partition, int sample_size,
    int min_per_partition, std::mt19937_64& gen) {
  std::vector<T> sampled;
  std::vector<T> unsampled;

  // First get min_per_partition from each partition.
  for (auto& partition_elements : partition) {
    int index =
        InPlaceReservoirSample(min_per_partition, gen, partition_elements);
    sampled.insert(sampled.end(), partition_elements.begin(),
                   partition_elements.begin() + index);
    unsampled.insert(unsampled.end(), partition_elements.begin() + index,
                     partition_elements.end());
  }

  // Complete the rest of the sample from the unsampled elements.
  int remaining_to_sample = sample_size - sampled.size();
  // If remaining_to_sample is negative, it means our threshold per allele
  // results in more than sample_size.
  if (remaining_to_sample < 0) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Threshold of ", min_per_partition,
        " per partition results in more than sample_size elements."));
  }
  int index = InPlaceReservoirSample(remaining_to_sample, gen, unsampled);
  sampled.insert(sampled.end(), unsampled.begin(), unsampled.begin() + index);
  return sampled;
}

}  // namespace learning::genomics::deepvariant::sampling

#endif  // LEARNING_GENOMICS_DEEPVARIANT_SAMPLING_UTIL_H_
