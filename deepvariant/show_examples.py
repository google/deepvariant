# Copyright 2020 Google LLC.
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
"""Generate human-readable images from DeepVariant example pileups.

# Only --examples is required.

show_examples
  --examples /path/to/make_examples.tfrecord@64.gz
  --vcf /path/to/variants_of_interest.vcf
  --regions "4:10-100 5:400-500" # or the file name(s) of BED/BEDPE files
  --output /path/to/output_prefix
  --image_type both
  --num_records 200
  --verbose
"""

import gzip
import json
import os
from typing import Any, Callable, Optional, Sequence, Set

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import ranges
from third_party.nucleus.util import vis

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'examples', None, 'Path to a make_examples tfrecord file or '
    'many sharded files using e.g. make_examples.tfrecord@64.gz. '
    'May be gzipped.')
flags.DEFINE_string(
    'example_info_json', None, 'Path to one *example_info.json file containing '
    'the information of the channels for the examples.')
flags.DEFINE_string(
    'vcf', None, '[optional] Path to vcf file to filter by. '
    'This will output exclusively the loci that match a row in '
    'the VCF file by chromosome, position, and reference bases. '
    'The VCF can be headerless, so for example, running grep on a hap.py '
    'output VCF file to get all false positives. '
    'The VCF may be gzipped or uncompressed.')
flags.DEFINE_string('output', None,
                    '[optional] Output prefix to write image files to.')
flags.DEFINE_enum('image_type', 'channels', ['channels', 'RGB', 'both'],
                  'Show examples as "channels", "RGB", or "both".')
flags.DEFINE_string(
    'regions', None, '[optional] Space-separated list of regions to filter to. '
    'Elements can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE '
    'files. Coordinates are 1-based, like in the VCF.')
flags.DEFINE_integer('num_records', -1,
                     'Maximum number of loci to output images for.')
flags.DEFINE_bool(
    'annotation', True, 'Label images with channel labels and mark midpoints. '
    'True by default. Use --noannotation to turn off.')
flags.DEFINE_bool('verbose', False,
                  'Show ID for each example as images are created.')
flags.DEFINE_bool(
    'truth_labels', True, 'For examples with truth labels, '
    'add the truth label to the file name.')
flags.DEFINE_string(
    'column_labels', None, 'Comma-separated column labels to print on image. '
    'Defaults to the standard channel names of DeepVariant. '
    'Use --noannotation to remove them entirely.')
flags.DEFINE_integer('scale', 1, 'Scale image outputs x times.')

UPDATE_EVERY_N_EXAMPLES = 10000
MAX_SIZE_TO_PRINT = 5


def parse_vcf(vcf_path: str) -> Set[str]:
  """Parse VCF to extract a dict keyed by locus IDs.

  Args:
      vcf_path: string, a path to a VCF file, that is gzipped (.gz) or not.

  Returns:
      All locus IDs from the VCF, where each locus ID has form "chr:start_end".
  """
  # Read gzipped file or uncompressed file.
  if vcf_path.endswith('.gz'):
    vcf_reader = gzip.open(vcf_path)
  else:
    vcf_reader = open(vcf_path, 'r')

  # This is not using nucleus.io.vcf VcfReader because it needs to support
  # pieces of vcf files without headers.
  ids_from_vcf = set()
  for l in vcf_reader:
    if isinstance(l, bytes):
      l = l.decode('utf-8')
    if not l.startswith('#'):
      cols = l.strip().split()
      if len(cols) < 4:
        continue
      # VCF uses 1-based positions while Nucleus' variant proto is 0-based.
      # Subtracting 1 here converts the VCF positions to 0-based coordinates.
      pos = int(cols[1]) - 1
      # Format: chrom:start_refBases, e.g. X:100_A.
      locus_id = '{}:{}_{}'.format(cols[0], pos, cols[3])
      ids_from_vcf.add(locus_id)
  vcf_reader.close()
  return ids_from_vcf


def get_full_id(variant: variants_pb2.Variant, indices: Sequence[int]) -> str:
  alt_genotype_string = '|'.join([variant.alternate_bases[i] for i in indices])
  return '{}:{}_{}->{}'.format(variant.reference_name, variant.start,
                               variant.reference_bases, alt_genotype_string)


def get_short_id(variant: variants_pb2.Variant, indices: Sequence[int]) -> str:
  """Prepare a locus ID, shortening ref and alt if necessary."""

  pos_prefix = '{}:{}'.format(variant.reference_name, variant.start)

  ref_bases = variant.reference_bases

  alts = variant.alternate_bases

  # If any ref or alt strings are too long, shorten them all.
  if len(ref_bases) > MAX_SIZE_TO_PRINT or any(
      [len(alts[i]) > MAX_SIZE_TO_PRINT for i in indices]):

    # If any alts are the same length (rare but possible), include their IDs.
    use_alt_indices = len(set([len(a) for a in alts])) < len(alts)
    alt_types = []
    for i in indices:
      a = alts[i]
      diff = len(a) - len(ref_bases)
      optional_id = 'alt{}'.format(i) if use_alt_indices else ''
      alt_type = ''
      if diff < 0:
        alt_type = 'DEL{}bp'.format(-1 * diff)
      elif diff > 0:
        alt_type = 'INS{}bp'.format(diff)
      else:
        alt_type = 'MNP{}bp'.format(len(a))
      alt_types.append('{}{}'.format(optional_id, alt_type))
    return '{}_{}'.format(pos_prefix, '|'.join(alt_types))
  else:
    # If ref and alts are short enough, show them fully: e.g. A->AG|C
    alt_strings = [alts[i] for i in indices]
    return '{}_{}->{}'.format(pos_prefix, ref_bases, '|'.join(alt_strings))


def get_label(example: tf.train.Example) -> Optional[int]:
  val = example.features.feature['label'].int64_list.value
  if val:
    return int(val[0])
  else:
    return None


def create_region_filter(region_flag_string: str,
                         verbose: bool = False) -> Callable[[Any], Any]:
  """Create a function that acts as a regions filter.

  Args:
    region_flag_string: string from --regions.
    verbose: bool. Whether to print regions after parsing.

  Returns:
    A function that given a variant will return True or False whether the
        variant falls inside the regions.

  """
  if isinstance(region_flag_string, str):
    region_args = region_flag_string.split()
  regions = ranges.RangeSet.from_regions(region_args)
  if verbose:
    logging.info('Regions to filter to: %s',
                 ', '.join([ranges.to_literal(r) for r in regions]))

  def passes_region_filter(variant):
    for r in regions:
      if ranges.position_overlaps(variant.reference_name, variant.start, r):
        return True
    return False

  return passes_region_filter


def run():
  """Create pileup images from examples, filtered in various ways."""
  with errors.clean_commandline_error_exit():
    if FLAGS.column_labels and FLAGS.example_info_json:
      raise ValueError(
          'Set at most one of --column_labels or --example_info_json.')

    if FLAGS.column_labels:
      column_labels = FLAGS.column_labels.split(',')
    else:
      column_labels = None

    if FLAGS.example_info_json:
      example_info = json.load(tf.io.gfile.GFile(FLAGS.example_info_json, 'r'))
      column_labels = [
          deepvariant_pb2.DeepVariantChannelEnum.Name(x)
          for x in example_info['channels']
      ]

    filter_to_vcf = FLAGS.vcf is not None
    if filter_to_vcf:
      ids_from_vcf = parse_vcf(FLAGS.vcf)
      logging.info(
          'Found %d loci in VCF. '
          'Only examples matching these loci will be output.',
          len(ids_from_vcf))

    filter_to_region = FLAGS.regions is not None
    if filter_to_region:
      passes_region_filter = create_region_filter(
          region_flag_string=FLAGS.regions, verbose=FLAGS.verbose)

    # Use nucleus.io.tfrecord to read all shards.
    dataset = tfrecord.read_tfrecords(FLAGS.examples)

    # Check flag here to avoid expensive string matching on every iteration.
    make_rgb = FLAGS.image_type in ['both', 'RGB']
    make_channels = FLAGS.image_type in ['both', 'channels']

    num_scanned = 0
    num_output = 0
    for example in dataset:
      num_scanned += 1
      # Only when scanning many examples, print a dot for each one to
      # indicate that the script is making progress and not stalled.
      if num_scanned % UPDATE_EVERY_N_EXAMPLES == 0:
        if num_scanned == UPDATE_EVERY_N_EXAMPLES:
          print('Reporting progress below. Writing one dot every time {} '
                'examples have been scanned:'.format(UPDATE_EVERY_N_EXAMPLES))
        # Print another dot on the same line, using print since logging does
        # not support printing without a newline.
        print('.', end='', flush=True)

      # Extract variant from example.
      variant = vis.variant_from_example(example)
      locus_id = vis.locus_id_from_variant(variant)
      indices = vis.alt_allele_indices_from_example(example)

      # Optionally filter to variants in the VCF.
      if filter_to_vcf:
        # Check if the locus is in the VCF.
        if locus_id not in ids_from_vcf:
          # Skip this example since it doesn't match the VCF.
          continue

      if filter_to_region and not passes_region_filter(variant):
        continue

      # Use locus ID in the filename, replacing long alleles with INS/DEL sizes.
      locus_with_alt_id = get_short_id(variant, indices)

      # Examples of long alleles replaced with their sizes:
      # 20:62456134_INS103bp.png
      # 20:62481177_DEL51bp.png

      # Examples of short alleles where the full string is included:
      # 1:55424995_TC->T.png
      # 1:55424996_CT->CTT.png
      # 1:55424996_CT->C.png
      # 1:55424996_CT->TTT.png
      # 1:55424996_CT->C|CTT.png

      if FLAGS.verbose:
        logging.info('\nOutputting image for: %s', locus_with_alt_id)
        full_id = get_full_id(variant, indices)
        if locus_with_alt_id != full_id:
          logging.info(
              'ID above was shortened due to long ref/alt strings. '
              'Original: %s', full_id)

      # If the example has a truth label, optionally include it.
      optional_truth_label = ''
      if FLAGS.truth_labels:
        truth_label = get_label(example)
        if truth_label is not None:
          optional_truth_label = '_label{}'.format(truth_label)

      # Extract and format example into channels.
      channels = vis.channels_from_example(example)
      if column_labels is not None and len(column_labels) != len(channels):
        raise ValueError(
            '--column_labels must have {} names separated by commas, since '
            'there are {} channels in the examples. '
            'However, {} column labels were found: {}'.format(
                len(channels), len(channels), len(column_labels),
                ','.join(['"{}"'.format(x) for x in column_labels])))

      output_prefix = '{}_'.format(
          FLAGS.output) if FLAGS.output is not None else ''

      # Create directory for images if does not exist.
      if output_prefix:
        output_dir = os.path.dirname(output_prefix)
        if not tf.io.gfile.exists(output_dir):
          tf.io.gfile.makedirs(output_dir)

      # Create image with a grey-scale row of channels and save to file.
      if make_channels:
        channels_output = '{}channels_{}{}.png'.format(output_prefix,
                                                       locus_with_alt_id,
                                                       optional_truth_label)

        vis.draw_deepvariant_pileup(
            channels=channels,
            path=channels_output,
            scale=FLAGS.scale,
            show=False,
            annotated=FLAGS.annotation,
            labels=column_labels)

      # Create RGB image and save to file.
      if make_rgb:
        rgb_output = '{}rgb_{}{}.png'.format(output_prefix, locus_with_alt_id,
                                             optional_truth_label)
        vis.draw_deepvariant_pileup(
            channels=channels,
            composite_type='RGB',
            path=rgb_output,
            scale=FLAGS.scale,
            show=False,
            annotated=FLAGS.annotation,
            labels=column_labels)

      # Check if --num_records quota has been hit yet.
      num_output += 1
      if FLAGS.num_records != -1 and num_output >= FLAGS.num_records:
        break

    logging.info('Scanned %d examples and output %d images.', num_scanned,
                 num_output)

    if num_scanned == 0 and FLAGS.examples.startswith('gs://'):
      if sharded_file_utils.is_sharded_file_spec(FLAGS.examples):
        paths = sharded_file_utils.generate_sharded_filenames(FLAGS.examples)
        special_gcs_message = ('WARNING: --examples sharded files are either '
                               'all empty or do not exist. Please check that '
                               'the paths are correct:\n')
        for p in paths[0:3]:
          special_gcs_message += 'gsutil ls {}\n'.format(p)
        logging.warning(special_gcs_message)
      else:
        logging.warning(
            'WARNING: --examples file is either empty or does not exist. '
            'Please check that the path is correct: \n'
            'gsutil ls %s', FLAGS.examples)


def main(argv):
  logging.set_stderrthreshold('info')
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: show_examples.py does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv[1:])), errors.CommandLineError)
    run()


if __name__ == '__main__':
  flags.mark_flags_as_required(['examples'])
  app.run(main)
