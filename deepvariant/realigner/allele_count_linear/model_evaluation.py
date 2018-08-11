# Copyright 2018 Google Inc.
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
"""Outputs precision-recall for a sklearn model using AlleleCount features.

A special script is needed because an INDEL can be realigned even if it is
not labelled as such, as long as a position close to it is realigned.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import csv
from absl import app
from absl import flags
import pandas as pd
import six.moves
from sklearn.externals import joblib
import tensorflow as tf

from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'truth_variants', None,
    'Required. Tabix-indexed VCF file containing the truth variant calls for '
    'this labels which we use to label our examples.')
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call.')
flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
flags.DEFINE_string('input_model_pckl', None,
                    'Required. Path to an sklearn model in pckl format.')
flags.DEFINE_string(
    'eval_region', '20', 'Required. Region to evaluate on, in '
    'in the "chr:start-end", "chr:position" or "chr" format.')
flags.DEFINE_string('output_report_csv', None,
                    'Required. File path to write the precision recall data.')

_WINDOW_SIZE = 10000
_ALLELE_OP_STR = {
    deepvariant_pb2.SUBSTITUTION: 'SUBSTITUTION',
    deepvariant_pb2.INSERTION: 'INSERTION',
    deepvariant_pb2.DELETION: 'DELETION',
    deepvariant_pb2.SOFT_CLIP: 'SOFT_CLIP',
}
_THRESHOLDS = [float(x) / 2 for x in range(-10, 20)]
_RECALL_WINDOW = 30


def allele_count_to_vector(allele_count):
  row = collections.defaultdict(float)
  for operation in _ALLELE_OP_STR.values():
    row[operation] = 0
  row['ref_nonconfident_read_count'] = allele_count.ref_nonconfident_read_count
  row['ref_supporting_read_count'] = allele_count.ref_supporting_read_count
  for _, allele in allele_count.read_alleles.iteritems():
    row[_ALLELE_OP_STR[allele.type]] += allele.count
  return row


def _check_allele_count_quality(allele_count):
  return (allele_count.ref_base in 'ACGT' and
          (len(allele_count.read_alleles) or
           allele_count.ref_nonconfident_read_count != 0 or
           allele_count.ref_supporting_read_count != 0))


def compute_effective_recall(model, true_indels, sam_reader, ref_reader,
                             allele_counter_options, thresholds):
  """Window size aware recall computation.

  During realignment, the region around selected positions will be realigned.
  This means that if a position is marked for realignment close enough to an
  unrecognized INDEL, that INDEL will still be realigned.
  Since scikit-learn and pandas do not offer a feature for that definition of
  recall, we need to manually check the neighborhood of each INDEL.

  Args:
    model: a scikit-learn model implementing the decision_function method.
    true_indels: a list of nucleus.Variants.
    sam_reader: a nucleus.io.SamReader.
    ref_reader: a nucleus.io.RefFastaReader.
    allele_counter_options: a deepvariant.AlleleCounterOptions.
    thresholds: a list of threshold to compute recall on.

  Returns:
    A dict with a threshold as key and recall at that threshold as value.
  """
  scores = collections.defaultdict(float)
  n_indels = len(true_indels)

  for indel in true_indels:
    region = ranges.make_range(indel.reference_name,
                               indel.start - _RECALL_WINDOW,
                               indel.start + _RECALL_WINDOW)
    allele_counter = allelecounter.AlleleCounter(ref_reader.c_reader, region,
                                                 allele_counter_options)

    reads = list(sam_reader.query(region))
    for read in reads:
      allele_counter.add(read)
    counts = allele_counter.counts()

    x_region = pd.DataFrame.from_records([
        allele_count_to_vector(count)
        for count in counts
        if _check_allele_count_quality(count)
    ])
    x_region.fillna(0)

    y_score = model.decision_function(x_region[[
        'ref_nonconfident_read_count', 'ref_supporting_read_count',
        'SUBSTITUTION', 'INSERTION', 'DELETION', 'SOFT_CLIP'
    ]])

    # This INDEL will be realigned if any position within the window is
    # marked for realignment.
    y_max = y_score.max()
    for threshold in thresholds:
      if y_max > threshold:
        scores[threshold] += 1

  for score in scores:
    scores[score] /= n_indels
  return scores


def compute_precision(model, true_indels, sam_reader, ref_reader,
                      allele_counter_options, thresholds, eval_region):
  """Simple precision computation for a given set of thresholds.

  This implementation computes the precision in a sliding window fashion in
  order to limit the memory consumption since it can done on a large
  `eval_region`.

  Args:
    model: a scikit-learn model implementing the decision_function method.
    true_indels: a list of nucleus.Variants expected to come from `eval_region`.
    sam_reader: a nucleus.io.SamReader.
    ref_reader: a nucleus.io.RefFastaReader.
    allele_counter_options: a deepvariant.AlleleCounterOptions.
    thresholds: a list of threshold to compute recall on.
    eval_region: a nucleus.v1.Range of the range to evaluate on.

  Returns:
    A dict with a threshold as key and precision at that threshold as value.
  """
  # They all come from `eval_region` so they can be identified by start.
  indels = set([var.start for var in true_indels])
  true_positives = collections.defaultdict(float)
  positives = collections.defaultdict(float)

  for position in six.moves.range(eval_region.start, eval_region.end,
                                  _WINDOW_SIZE):
    region = ranges.make_range(eval_region.reference_name, position,
                               min(position + _WINDOW_SIZE, eval_region.end))
    allele_counter = allelecounter.AlleleCounter(ref_reader.c_reader, region,
                                                 allele_counter_options)

    reads = list(sam_reader.query(region))
    for read in reads:
      allele_counter.add(read)
    counts = allele_counter.counts()

    x_region_list = []
    for pos, count in enumerate(counts, start=position):
      if not _check_allele_count_quality(count):
        continue

      row = allele_count_to_vector(count)
      row['label'] = 1 if pos in indels else 0
      x_region_list.append(row)

    if not x_region_list:
      continue

    x_region = pd.DataFrame.from_records(x_region_list)
    x_region.fillna(0)

    y_score = model.decision_function(x_region[[
        'ref_nonconfident_read_count', 'ref_supporting_read_count',
        'SUBSTITUTION', 'INSERTION', 'DELETION', 'SOFT_CLIP'
    ]])

    positives_mask = x_region['label'] == 1
    for threshold in thresholds:
      positives[threshold] += sum(y_score > threshold)
      true_positives[threshold] += sum((y_score > threshold) & positives_mask)

    if position % 100000 == 0:
      print('processed %d positions out of %d -- %2f complete' %
            (position - eval_region.start, eval_region.end - eval_region.start,
             float(position - eval_region.start) /
             (eval_region.end - eval_region.start)))

  precisions = {}
  for threshold in thresholds:
    if positives[threshold] == 0:
      # We arbitrarily decide that if we find no positives our precision is NaN.
      precisions[threshold] = float('nan')
    else:
      precisions[threshold] = true_positives[threshold] / positives[threshold]
  return precisions


def model_evaluation_runner(truth_variants, reads, ref, input_model_pckl,
                            eval_region, output_report_csv):
  """Outputs precision-recall for a sklearn model using AlleleCount features.

  Args:
    truth_variants: path to the VCF.
    reads: path to the reads BAM.
    ref: path to the reference FASTA.
    input_model_pckl: path to read the LogisticRegression pickle from.
    eval_region: str, region to evaluate on in the 'chr:start-end',
      'chr:position' or 'chr' format.
    output_report_csv: path to the output report csv.

  Raises:
    ValueError: if eval_region cannot be parsed.
  """
  sam_reader = sam.SamReader(reads)
  ref_reader = fasta.RefFastaReader(ref)

  read_reqs = reads_pb2.ReadRequirements(
      min_base_quality=10,
      min_mapping_quality=10,
      min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT)
  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=1, read_requirements=read_reqs)

  model = joblib.load(input_model_pckl)

  with vcf.VcfReader(truth_variants) as vcf_reader:
    region = ranges.parse_literal(
        eval_region, contig_map=ranges.contigs_dict(ref_reader.header.contigs))
    true_indels = [
        var for var in vcf_reader.query(region) if (variant_utils.is_indel(var))
    ]

  precisions = compute_precision(model, true_indels, sam_reader, ref_reader,
                                 allele_counter_options, _THRESHOLDS, region)
  recalls = compute_effective_recall(model, true_indels, sam_reader, ref_reader,
                                     allele_counter_options, _THRESHOLDS)

  with tf.gfile.GFile(output_report_csv, 'w') as csvfile:
    fieldnames = ['threshold', 'precision', 'recall']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for threshold in _THRESHOLDS:
      writer.writerow({
          'threshold': threshold,
          'precision': precisions[threshold],
          'recall': recalls[threshold]
      })


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  model_evaluation_runner(
      truth_variants=FLAGS.truth_variants,
      reads=FLAGS.reads,
      ref=FLAGS.ref,
      input_model_pckl=FLAGS.input_model_pckl,
      eval_region=FLAGS.eval_region,
      output_report_csv=FLAGS.output_report_csv)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'ref', 'reads', 'truth_variants', 'input_model_pckl', 'output_report_csv'
  ])
  app.run(main)
