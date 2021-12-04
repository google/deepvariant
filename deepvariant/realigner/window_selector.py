# Copyright 2017 Google LLC.
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
"""Determine genomic ranges to perform local assembly."""

from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from deepvariant.protos import deepvariant_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner.python import window_selector as cpp_window_selector


def _candidates_from_reads(config, ref_reader, reads, region):
  """Returns a list of candidate positions.

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref_reader: GenomeReference. Indexed reference genome to query bases.
    reads: list[nucleus.protos.Read]. The reads we are processing into candidate
      positions.
    region: nucleus.protos.Range. The region we are processing.

  Returns:
    A list. The elements are reference positions within region.

  Raises:
    ValueError: if config.window_selector_model.model_type isn't a valid enum
    name in realigner_pb2.WindowSelectorModel.ModelType.
  """
  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      read_requirements=reads_pb2.ReadRequirements(
          min_mapping_quality=config.min_mapq,
          min_base_quality=config.min_base_quality),
      keep_legacy_behavior=config.keep_legacy_behavior)
  expanded_region = ranges.expand(
      region,
      config.region_expansion_in_bp,
      contig_map=ranges.contigs_dict(ref_reader.header.contigs))

  allele_counter = allelecounter.AlleleCounter(ref_reader.c_reader,
                                               expanded_region, [],
                                               allele_counter_options)

  for read in reads:
    allele_counter.add(read, 'placeholder_sample_id')

  model_type = config.window_selector_model.model_type
  if model_type == realigner_pb2.WindowSelectorModel.VARIANT_READS:
    return _variant_reads_threshold_selector(
        allele_counter, config.window_selector_model.variant_reads_model,
        expanded_region)
  elif model_type == realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR:
    return _allele_count_linear_selector(
        allele_counter, config.window_selector_model.allele_count_linear_model,
        expanded_region)
  else:
    raise ValueError('Unknown enum option "{}" for '
                     'WindowSelectorModel.model_type'.format(
                         config.window_selector_model.model_type))


def _variant_reads_threshold_selector(allele_counter, model_conf,
                                      expanded_region):
  """Returns a list of candidate positions.

  Following cigar operations generate candidate position:
    - ALIGNMENT_MATCH, SEQUENCE_MISMATCH, SEQUENCE_MATCH: at mismatch positions
      in the read when compared to the reference sequence.
    - DELETE: at positions within [cigar_start, cigar_start + cigar_len)
    - INSERT, CLIP_SOFT: at positions within
        [cigar_start - cigar_len, cigar_start + cigar_len)

   Note. Function implementation has changed to return positions beyond input
   region in case we have variants there. See the change at internal and
   internal.

  Args:
    allele_counter: learning.genomics.deepvariant.realigner.AlleleCounter in the
      considered region.
    model_conf: learning.genomics.deepvariant.realigner
      .WindowSelectorOptions.VariantReadsThresholdModel options determining the
      behavior of this window selector.
    expanded_region: nucleus.protos.Range. The region we are processing.

  Returns:
    A list. The elements are reference positions within region.
  """

  counts_vec = cpp_window_selector.variant_reads_candidates_from_allele_counter(
      allele_counter)

  return [
      expanded_region.start + i
      for i, count in enumerate(counts_vec)
      if (count >= model_conf.min_num_supporting_reads and
          count <= model_conf.max_num_supporting_reads)
  ]


def _allele_count_linear_selector(allele_counter, model_conf, expanded_region):
  """Returns a list of candidate positions.

  Candidate positions for realignment are generated by scoring each location.
  The score at a location is a weighted sum of the number of reads with each
  CIGAR operation at the location, where the weights are determined by the model
  coefficients. Locations whose score exceed the model decision boundary value
  are used to create realignment windows.

  Note. Function implementation has changed to return positions beyond input
  region in case we have variants there. See the change at internal and
  internal.

  Args:
    allele_counter: learning.genomics.deepvariant.realigner.AlleleCounter in the
      considered region.
    model_conf: learning.genomics.deepvariant.realigner
      .WindowSelectorOptions.AlleleCountLinearModel options determining the
      behavior of this window selector.
    expanded_region: nucleus.protos.Range. The region we are processing.

  Returns:
    A list. The elements are reference positions within region.
  """

  scores_vec = (
      cpp_window_selector.allele_count_linear_candidates_from_allele_counter(
          allele_counter, model_conf))

  return [
      expanded_region.start + i
      for i, score in enumerate(scores_vec)
      if score > model_conf.decision_boundary
  ]


def _candidates_to_windows(config, candidate_pos, ref_name):
  """"Process candidate positions to determine windows for local assembly.

  Windows are within range of
    [min(pos) - config.min_windows_distance,
     max(pos) + config.min_windows_distance)

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    candidate_pos: A list of ref_pos.
    ref_name: Reference name, used in setting the output
      genomics.range.reference_name value.

  Returns:
    A sorted list of nucleus.protos.Range protos for all windows in this region.
  """
  windows = []

  def _add_window(start_pos, end_pos):
    windows.append(
        ranges.make_range(ref_name, start_pos - config.min_windows_distance,
                          end_pos + config.min_windows_distance))

  start_pos, end_pos = None, None
  for pos in sorted(candidate_pos):
    if start_pos is None:
      start_pos = pos
      end_pos = pos
    # We need to check if the previous end_pos is within 2*window_distance as we
    # generate a window of radius window_distance around each position.
    #
    #   <-------end_pos------->
    #                          <-------pos------->
    # where window_distance = ------->
    #
    # If this is the case, we need to merge the two windows.
    elif pos > end_pos + 2 * config.min_windows_distance:
      _add_window(start_pos, end_pos)
      start_pos = pos
      end_pos = pos
    else:
      end_pos = pos
  if start_pos is not None:
    _add_window(start_pos, end_pos)

  return sorted(windows, key=ranges.as_tuple)


def select_windows(config, ref_reader, reads, region):
  """"Process reads to determine candidate windows for local assembly.

  Windows are within range of
    [0 - config.min_windows_distance, ref_len + config.min_windows_distance)

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref_reader: GenomeReference. Indexed reference genome to query bases.
    reads: A list of genomics.Read records.
    region: nucleus.protos.Range. The region we are processing.

  Returns:
    A list of nucleus.protos.Range protos sorted by their genomic position.
  """
  # This is a fast path for the case where we have no reads, so we have no
  # windows to assemble.
  if not reads:
    return []

  candidates = _candidates_from_reads(config, ref_reader, reads, region)
  return _candidates_to_windows(config, candidates, region.reference_name)
