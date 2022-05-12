# Copyright 2018 Google LLC.
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

# pylint: disable=line-too-long
r"""Shuffle tf.Example files using beam.

To run locally:
1) Install beam on your machine following the instructions at
   https://beam.apache.org/get-started/quickstart-py/

2) Copy any inputs to be on local disk.

3) Run
  python path/to/shuffle_tfrecords_beam.py \
    --input_pattern_list="/tmp/some.examples-?????-of-00200.tfrecord.gz" \
    --output_pattern_prefix="/tmp/training.examples" \
    --output_dataset_name="HG001" \
    --runner=DirectRunner

To run on Google Cloud Dataflow Service:
1) Follow the Google Cloud Dataflow setup instructions at
https://beam.apache.org/documentation/runners/dataflow/

2) Upload this file to your GCE instance.

3) Run
  python shuffle_tfrecords_beam.py \
  --job_name=shuffle-tfrecords \
  --input_pattern_list="gs://YOUR_INPUT_BUCKET/A.tfrecord.gz" \
  --output_pattern_prefix="gs://YOUR_OUTPUT_BUCKET/training.examples" \
  --output_dataset_name="HG001" \
  --runner=DataflowRunner \
  --project=SET_YOUR_PROJECT_ID_HERE \
  --staging_location=gs://YOUR_BUCKET_NAME/AND_STAGING_DIRECTORY \
  --temp_location=gs://YOUR_BUCKET_NAME/AND_TEMP_DIRECTORY

4) (Optional) To monitor or cancel the job while it is running, you can
use either the Dataflow Monitoring Interface
https://cloud.google.com/dataflow/pipelines/dataflow-monitoring-intf
or the Dataflow Command-line Interface
https://cloud.google.com/dataflow/pipelines/dataflow-command-line-intf
"""
# pylint: enable=line-too-long

import argparse
import hashlib
import logging
import os
import textwrap

import apache_beam as beam
from apache_beam import coders
from apache_beam.options.pipeline_options import PipelineOptions
import tensorflow.compat.v1 as tf

COMMENT_HEADER = """#
# --input_pattern_list={}
# --output_pattern_prefix={}
#
"""


def parse_cmdline(argv):
  """Parse the commandline into known and pipeline arguments.

  The known arguments are required for this specific program to function,
  and the other pipeline arguments can be used to configure beam and the
  specific beam backend being used.  See
  https://github.com/apache/beam/blob/master/sdks/python/apache_beam/options/pipeline_options.py
  for a list and description of the pipeline arguments accepted.

  Args:
    argv: List containing command-line arguments.

  Returns:
    A pair, the first of which are the known (non-pipeline) arguments
    and the second of which are the pipeline arguments.
  """
  parser = argparse.ArgumentParser()

  parser.add_argument(
      '--input_pattern_list',
      help='Comma-separated list of TFRecord filename patterns.')
  parser.add_argument(
      '--output_pattern_prefix',
      help='Filename pattern for the output TFRecords.')
  parser.add_argument(
      '--output_dataset_config_pbtxt',
      help='Optional.  If set, print out a human-readable version of '
      'DeepVariantDatasetConfig.')
  parser.add_argument(
      '--output_dataset_name',
      help='Optional unless --output_dataset_config_pbtxt is set.')

  known_args, pipeline_args = parser.parse_known_args(argv)

  return known_args, pipeline_args


def read_from_tfrecords_files(pipeline, input_filename_pattern_list):
  """Reads records from TFRecord files.

  Args:
    pipeline: Beam pipeline object.
    input_filename_pattern_list: List of filename patterns.

  Returns:
    A PCollection of read tf.Examples.
  """
  readers = []
  for i, filepattern in enumerate(input_filename_pattern_list):
    readers.append(pipeline
                   | 'ReadTFRecordFiles_{}[{}]'.format(i, filepattern) >> beam
                   .io.ReadFromTFRecord(filepattern, coder=coders.BytesCoder()))
  return readers | 'Flatten' >> beam.Flatten()


def shuffle_records(input_examples):
  """Shuffles the input_examples in a effectively random order."""

  def sha1(input_bytes):
    """Returns the sha1 hash of input_bytes."""
    m = hashlib.sha1()
    m.update(input_bytes)
    return m.digest()

  return (input_examples
          | 'Randomize' >> beam.Map(lambda x: (sha1(x), x))
          | 'Groupby' >> beam.GroupByKey()
          | 'DropKey' >> beam.FlatMap(lambda x: x[1]))


def count_records_per_label(input_examples):
  """Shuffles the input_examples in a effectively random order."""

  def label_example(input_bytes):
    """Returns the label of input_example."""
    example = tf.train.Example.FromString(input_bytes)
    label = example.features.feature['label'].int64_list.value[0]
    return label

  return (
      input_examples
      | 'LabelExample' >> beam.Map(lambda x: (label_example(x), x))
      | 'CountPerLabel' >> beam.combiners.Count.PerKey()
      |
      'ToString' >> beam.Map(lambda kv: u'# class{}: {}\n'.format(kv[0], kv[1]))
      | 'Concat1' >> beam.CombineGlobally(''.join))


def make_config_string(name, tfrecord_path, num_examples):
  return textwrap.dedent("""
  name: "{}"
  tfrecord_path: "{}-?????-of-?????.tfrecord.gz"
  num_examples: {}
  """.format(name, tfrecord_path, num_examples))


def write_summary_string_to_file(pipeline, output_examples, input_pattern_list,
                                 dataset_name, output_pattern_prefix,
                                 output_filename):
  """Writes a file summarizing the PCollection of Examples.

  Args:
    pipeline: Beam pipeline object.
    output_examples: PCollection of examples.
    input_pattern_list: str. A comma-separated string of input files.
    dataset_name: str. The name of the dataset to be written in the output.
    output_pattern_prefix: str. The prefix of the sharded output files.
    output_filename: the output text file that contains the summary that can be
      parsed into DeepVariantDatasetConfig.
  """

  # Beam currently has no way to materialize pipeline values, so we have
  # to construct the file entirely in Beam pipeline operations.
  comment_str = pipeline | 'CreateFileHeader' >> beam.Create(
      [COMMENT_HEADER.format(input_pattern_list, output_pattern_prefix)])
  num_examples = (
      output_examples
      | 'CountOutputExamples' >> beam.combiners.Count.Globally())
  config_str = num_examples | 'MakeConfigStr' >> beam.Map(
      lambda n: make_config_string(dataset_name, output_pattern_prefix, n))

  num_examples_by_labels = count_records_per_label(output_examples)
  merged_strings = (comment_str, num_examples_by_labels,
                    config_str) | 'FlattenStrs' >> beam.Flatten()
  _ = (
      merged_strings
      | 'Concat2' >> beam.CombineGlobally(''.join)
      | 'WriteToFile' >> beam.io.WriteToText(
          output_filename,
          shard_name_template='',
          header='# Generated by shuffle_tfrecords_beam.py'))


def main(argv=None):
  """Main entry point; defines and runs the pipeline."""
  known_args, pipeline_args = parse_cmdline(argv)
  # Copy over the example_info.json file before the pipeline starts.
  example_info_json = '{}*example_info.json'.format(
      known_args.input_pattern_list.split(',')[0])
  example_info_json_list = tf.io.gfile.glob(example_info_json)
  if example_info_json_list and tf.io.gfile.exists(example_info_json_list[0]):
    training_dir = os.path.dirname(known_args.output_pattern_prefix)
    if not tf.io.gfile.exists(training_dir):
      tf.io.gfile.makedirs(training_dir)
    output_example_info_json = os.path.join(training_dir, 'example_info.json')
    if not tf.io.gfile.exists(output_example_info_json):
      tf.io.gfile.copy(example_info_json_list[0], output_example_info_json)

  pipeline_options = PipelineOptions(pipeline_args)
  with beam.Pipeline(options=pipeline_options) as p:
    input_examples = read_from_tfrecords_files(
        p, known_args.input_pattern_list.split(','))
    output_examples = shuffle_records(input_examples)

    _ = output_examples | beam.io.WriteToTFRecord(
        file_path_prefix=known_args.output_pattern_prefix,
        file_name_suffix='.tfrecord.gz',
        coder=coders.BytesCoder())
    if known_args.output_dataset_config_pbtxt:
      if not known_args.output_dataset_name:
        raise ValueError('Need to set output_dataset_name.')
      write_summary_string_to_file(p, output_examples,
                                   known_args.input_pattern_list,
                                   known_args.output_dataset_name,
                                   known_args.output_pattern_prefix,
                                   known_args.output_dataset_config_pbtxt)


if __name__ == '__main__':
  logging.getLogger().setLevel(logging.INFO)
  main()
