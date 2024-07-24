# Copyright 2021 Google LLC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Module for partitioning contigs."""

from typing import Optional, Sequence, Union



from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.util import ranges


def parse_regions_flag(
    regions_flag_value: Union[str, Sequence[str]],
) -> Sequence[str]:
  """Splits a regions flag value into a list of regions if necessary."""
  if isinstance(regions_flag_value, str):
    regions_flag_value = regions_flag_value.split()
  return regions_flag_value


def build_calling_regions(
    contigs: Sequence[reference_pb2.ContigInfo],
    regions_to_include: Sequence[str],
    regions_to_exclude: Sequence[str],
    ref_n_regions: Optional[Sequence[range_pb2.Range]],
) -> ranges.RangeSet:
  """Builds a RangeSet containing the regions we should call variants in.

  This function intersects the Ranges spanning all of the contigs with those
  from regions_to_include, if not empty, and removes all of the regions in
  regions_to_exclude.

  Args:
    contigs: Sequence of ContigInfo protos. Used to determine the initial ranges
      to process (i.e., all bases of these contigs).
    regions_to_include: RangeSet or iterable that can be converted to a
      RangeSet.
    regions_to_exclude: RangeSet or iterable that can be converted to a
      RangeSet.
    ref_n_regions: List of Range containing non DNA bases to exclude.

  Returns:
    A RangeSet.
  """
  # Initially we are going to call everything in the reference.
  regions = ranges.RangeSet.from_contigs(contigs)

  # If we provided a regions to include, intersect it with all of the regions,
  # producing a common set of regions between the reference and the provided
  # calling regions.
  contig_dict = ranges.contigs_dict(contigs)
  if regions_to_include:
    regions = regions.intersection(
        ranges.RangeSet.from_regions(regions_to_include, contig_dict)
    )

  if ref_n_regions:
    regions.exclude_regions(ranges.RangeSet(ref_n_regions))

  # If we provided regions to exclude, intersect those with the existing calling
  # regions to further refine our set of contigs to process.
  if regions_to_exclude:
    # exclude_regions mutates regions.
    regions.exclude_regions(
        ranges.RangeSet.from_regions(regions_to_exclude, contig_dict)
    )

  return regions


def partition_calling_regions(
    calling_regions: ranges.RangeSet, num_partitions: int
) -> Sequence[Sequence[range_pb2.Range]]:
  """Splits the calling regions into N contiguous partitions.

  This is a simpler version of https://en.wikipedia.org/wiki/Partition_problem,
  since we have to preserve order. This algorithm runs in 3 stages:
    1. Calculate the partition size across the entire genome by BP count.
    2. Group the partitions greedily such that all groups exceed the partition
       size. This will produce *fewer* groups than `num_partitions`.
    3. Keep splitting the largest groups in half until the target number of
       partition is reached.
    4. Sort everything such that we preserve the order in `contigs`.

  Args:
    calling_regions: the set of all calling regions to partition.
    num_partitions: number of partitions to produce.

  Returns:
    partition_groups: list of groups of nucleus.genomics.v1.Range.
  """
  # Split into target partition size.
  total_bps = sum(r.end - r.start for r in calling_regions)
  max_partition_size = total_bps // num_partitions
  partitions = list(calling_regions.partition(max_partition_size))

  # Group the partitions such that no group exceeds the max size.
  partition_groups = []
  current_group = []
  for partition in partitions:
    if sum(p.end - p.start for p in current_group) > max_partition_size:
      partition_groups.append(current_group)
      current_group = []
    current_group.append(partition)
  if current_group:
    partition_groups.append(current_group)

  # The above will produce fewer groups than `num_partitions`, because
  # all groups grew past the target size. We fix this now by splitting
  # the largest groups in half until we have `num_partitions` groups.
  while len(partition_groups) < num_partitions:
    partition_groups.sort(key=lambda ps: sum(p.end - p.start for p in ps))
    largest_partition = partition_groups.pop()
    mid_point = len(largest_partition) // 2
    ps_1 = largest_partition[:mid_point]
    ps_2 = largest_partition[mid_point:]
    partition_groups.extend([ps_1, ps_2])

  # Sort the partitions into their initial config
  partition_groups.sort(key=lambda ps: partitions.index(ps[0]))
  return partition_groups


def build_and_partition_calling_regions(
    contigs: Sequence[reference_pb2.ContigInfo],
    num_partitions: int,
    calling_regions: str,
) -> Sequence[Sequence[range_pb2.Range]]:
  """Builds and partitions the calling regions."""
  regions = build_calling_regions(
      contigs=contigs,
      regions_to_include=parse_regions_flag(calling_regions),
      regions_to_exclude=[],
      ref_n_regions=[],
  )
  return partition_calling_regions(regions, num_partitions)
