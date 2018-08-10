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
"""Generates a trained AlleleCountLinearModel from given VCF and BAM."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import random
from absl import app
from absl import flags
import pandas as pd
from sklearn.externals import joblib
import sklearn.linear_model
import sklearn.metrics
from sklearn.utils import shuffle
import tensorflow as tf

from deepvariant.protos import deepvariant_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.python import allelecounter
from google.protobuf import text_format
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
    'this sample which we use to label our examples.')
flags.DEFINE_string('from_contig', '1',
                    'Contig from which the baseline will be sampled.')
flags.DEFINE_string(
    'exclude_contig', '20',
    'Contig to be excluded from the training so that it can be used for '
    'evaluation.')
flags.DEFINE_string(
    'output_model_proto', None,
    'Path to write the AlleleCountLinearModel proto in text format.')
flags.DEFINE_string(
    'output_model_pckl', None,
    'Path to write the trained LogisticRegression in pickle format.')
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call.')
flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
flags.DEFINE_integer('random_seed', 42, 'Random seed used for reproducibility.')
flags.DEFINE_float(
    'indel_weight', 3,
    'How much to weight indels relative to other classes in training.')

_REF_LABEL = 0
_SNP_LABEL = 0
_INDEL_LABEL = 1
_ALLELE_OP_STR = {
    deepvariant_pb2.UNSPECIFIED: 'UNSPECIFIED',
    deepvariant_pb2.REFERENCE: 'REFERENCE',
    deepvariant_pb2.SUBSTITUTION: 'SUBSTITUTION',
    deepvariant_pb2.INSERTION: 'INSERTION',
    deepvariant_pb2.DELETION: 'DELETION',
    deepvariant_pb2.SOFT_CLIP: 'SOFT_CLIP',
}
# Empirically found to be the best parameter on HG0001.
_DEFAULT_THRESHOLD = 3

PositionWrapper = collections.namedtuple('PositionWrapper',
                                         ['reference_name', 'start', 'label'])


def _is_valid(allele_count):
  return (allele_count.ref_base in 'ACGT' and
          (len(allele_count.read_alleles) or
           allele_count.ref_nonconfident_read_count != 0 or
           allele_count.ref_supporting_read_count != 0))


def _position_to_features(sam_reader, allele_counter, region, position,
                          exclude_contig):
  """Extracts the AlleleCount data from a given position."""
  # We build the AlleleCount at the position.
  reads = sam_reader.query(region)
  for read in reads:
    allele_counter.add(read)
  counts = allele_counter.counts()
  assert len(counts) == 1
  allele_count = counts[0]

  if not _is_valid(allele_count) or (position.reference_name == exclude_contig):
    return None

  # We turn that AlleleCount into a vector feedable to scikitlearn.
  row = dict()
  row['ref_nonconfident_read_count'] = allele_count.ref_nonconfident_read_count
  row['ref_supporting_read_count'] = allele_count.ref_supporting_read_count
  # We need to make sure that all columns exist.
  for operation in _ALLELE_OP_STR.values():
    row[operation] = 0
  for _, allele in allele_count.read_alleles.iteritems():
    row[_ALLELE_OP_STR[allele.type]] += allele.count
  row['label'] = position.label
  row['reference_name'] = position.reference_name
  row['position'] = position.start

  return row


def generate_positions(vcf_reader, ref_reader, baseline_contig):
  """Gets all INDELs position and an equal amount of SNPs and random positions.

  Args:
    vcf_reader: a nucleus.io.VcfReader.
    ref_reader: a nucleus.io.RefFastaReader.
    baseline_contig: contig from which to sample baseline positions.

  Returns:
    A list of PositionWrapper.
  """
  variants = [variant for variant in vcf_reader]
  indels_positions = [
      PositionWrapper(var.reference_name, var.start, _INDEL_LABEL)
      for var in variants
      if variant_utils.is_indel(var)
  ]
  n_indels = len(indels_positions)

  # We sort by position for better data locality.
  snps = [var for var in variants if variant_utils.is_snp(var)]
  snps_positions = [
      PositionWrapper(var.reference_name, var.start, _SNP_LABEL)
      for var in random.sample(snps, min(len(snps), n_indels))
  ]

  contig_size = ref_reader.contig(baseline_contig).n_bases
  # NOTE: Though unlikely, these random positions can end up on actual
  # variants.
  baseline_positions = [
      PositionWrapper(baseline_contig, pos, _REF_LABEL)
      for pos in random.sample(xrange(contig_size), min(contig_size, n_indels))
  ]

  return sorted(indels_positions + snps_positions + baseline_positions)


def generate_data(vcf_reader, ref_reader, sam_reader, baseline_contig,
                  exclude_contig):
  """Generates a pandas.DataFrame summarizing the AlleleCount at each position.

  The features included are:
        - 'ref_nonconfident_read_count'
        - 'ref_supporting_read_count'
        - 'SUBSTITUTION'
        - 'INSERTION'
        - 'DELETION'
        - 'SOFT_CLIP'
        - 'label'
  These features are extracted from the AlleleCount proto at the concerned
  position.

  Args:
    vcf_reader: a nucleus.io.VcfReader.
    ref_reader: a nucleus.io.RefFastaReader.
    sam_reader: a nucleus.io.SamReader.
    baseline_contig: string, contig from which to sample baseline positions.
    exclude_contig: string, contig to exclude for test purposes.

  Returns:
    pandas.Dataframe object.
  """

  # These parameters are the ones used in make_examples.
  read_reqs = reads_pb2.ReadRequirements(
      min_base_quality=10,
      min_mapping_quality=10,
      min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT)
  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=1, read_requirements=read_reqs)

  training_positions = generate_positions(vcf_reader, ref_reader,
                                          baseline_contig)
  positions_records = []

  for position in training_positions:
    region = ranges.make_range(position.reference_name, position.start,
                               position.start + 1)
    allele_counter = allelecounter.AlleleCounter(ref_reader.c_reader, region,
                                                 allele_counter_options)
    row = _position_to_features(sam_reader, allele_counter, region, position,
                                exclude_contig)
    if row is not None:
      positions_records.append(row)

  df = pd.DataFrame(positions_records)
  df = df.fillna(0)
  df = shuffle(df)
  return df


def train_model(data, indel_weight):
  """Generates a trained LogisticRegression model from an input dataframe."""
  x_train = data[[
      'ref_nonconfident_read_count', 'ref_supporting_read_count',
      'SUBSTITUTION', 'INSERTION', 'DELETION', 'SOFT_CLIP'
  ]]
  y_train = data['label']
  model = sklearn.linear_model.LogisticRegression(
      tol=1e-7,
      class_weight={
          _INDEL_LABEL: indel_weight,
          _REF_LABEL: 1
      },
      penalty='l2')
  model.fit(x_train, y_train)
  return model


def model_to_proto(model):
  """Returns an allele count-based linear WindowSelectorModel."""
  # AFAIK sklearn does not provide a way to extract the coefficients based on
  # the columns of its input, we thus have to use the fact we know their order.
  allele_count_linear_model = (
      realigner_pb2.WindowSelectorModel.AlleleCountLinearModel(
          bias=model.intercept_[0],
          coeff_soft_clip=model.coef_[0][5],
          coeff_substitution=model.coef_[0][2],
          coeff_insertion=model.coef_[0][3],
          coeff_deletion=model.coef_[0][4],
          coeff_reference=model.coef_[0][1],
          decision_boundary=_DEFAULT_THRESHOLD))
  return realigner_pb2.WindowSelectorModel(
      model_type=realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR,
      allele_count_linear_model=allele_count_linear_model)


def generate_trained_model_runner(
    truth_variants, reads, ref, output_model_proto, output_model_pckl,
    exclude_contig, from_contig, random_seed, indel_weight):
  """Runner for generate_trained_model.

  Args:
    truth_variants: path to the VCF.
    reads: path to the reads BAM.
    ref: path to the reference FASTA.
    output_model_proto: path to write the AlleleCountLinearModel proto.
    output_model_pckl: path to write the LogisticRegression pickle.
    exclude_contig: string identifier of a contig to exclude from training,
    from_contig: string identifier of the contig from which we sample baseline.
    random_seed: int used as random seed for reproducibility.
    indel_weight: float of the weight od indels relative to the rest in
      the training.
  """
  vcf_reader = vcf.VcfReader(truth_variants)
  ref_reader = fasta.RefFastaReader(ref)
  sam_reader = sam.SamReader(reads)

  random.seed(random_seed)

  dataframe = generate_data(vcf_reader, ref_reader, sam_reader, from_contig,
                            exclude_contig)
  model = train_model(dataframe, indel_weight=indel_weight)

  if output_model_pckl:
    joblib.dump(model, output_model_pckl)

  model_proto = model_to_proto(model)
  with tf.gfile.GFile(output_model_proto, 'w') as f:
    f.write(text_format.MessageToString(model_proto))


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  generate_trained_model_runner(
      truth_variants=FLAGS.truth_variants,
      reads=FLAGS.reads,
      ref=FLAGS.ref,
      output_model_proto=FLAGS.output_model_proto,
      output_model_pckl=FLAGS.output_model_pckl,
      exclude_contig=FLAGS.exclude_contig,
      from_contig=FLAGS.from_contig,
      random_seed=FLAGS.random_seed,
      indel_weight=FLAGS.indel_weight)


if __name__ == '__main__':
  flags.mark_flags_as_required(
      ['truth_variants', 'ref', 'reads', 'output_model_proto'])
  app.run(main)
