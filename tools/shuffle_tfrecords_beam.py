# Copyright 2018 Google Inc.  All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License

# pylint: disable=line-too-long
r"""Shuffle tf.Example files using beam.

To run locally:
1) Install beam on your machine following the instructions at
   https://beam.apache.org/get-started/quickstart-py/

2) Copy any inputs to be on local disk.

3) Run
  python path/to/shuffle_tfrecords_beam.py \
    --input_pattern_list="/tmp/some.examples-?????-of-00200.tfrecord.gz" \
    --output_pattern="/tmp/training.examples.tfrecord.gz" \
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
  --output_pattern="gs://YOUR_OUTPUT_BUCKET/training.examples.tfrecord.gz" \
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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import hashlib
import logging

import apache_beam as beam
from apache_beam import coders
from apache_beam.options.pipeline_options import PipelineOptions


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
      '--output_pattern', help='Filename pattern for the output TFRecords.')
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
                   | 'ReadTFRecordFiles_{}[{}]'.format(i, filepattern) >> beam.
                   io.ReadFromTFRecord(filepattern, coder=coders.BytesCoder()))
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


def main(argv=None):
  """Main entry point; defines and runs the pipeline."""
  known_args, pipeline_args = parse_cmdline(argv)
  pipeline_options = PipelineOptions(pipeline_args)
  with beam.Pipeline(options=pipeline_options) as p:
    input_examples = read_from_tfrecords_files(
        p, known_args.input_pattern_list.split(','))

    output_examples = shuffle_records(input_examples)

    # redacted

    # WriteToTFRecord's default compression_type is "AUTO", which uses
    # the file path's extension to determine the compression type.
    # redacted
    # and shard_name_template.
    _ = output_examples | beam.io.WriteToTFRecord(
        known_args.output_pattern, coder=coders.BytesCoder())


if __name__ == '__main__':
  logging.getLogger().setLevel(logging.INFO)
  main()
