# pylint: disable=line-too-long
r"""Globally shuffle DeepVariant tf.Example TFRecords with plain PySpark.

This is a drop-in alternative to shuffle_tfrecords_beam.py for AWS EMR, any
other Spark cluster, or a single machine. It does not use Apache Beam, so there
is no Beam job server, no Docker SDK harness and no Maven download at runtime.

Memory use is bounded: input files are streamed record by record, and the
global shuffle is a Spark sort-based shuffle that spills to local disk. Neither
the driver nor any executor ever holds the whole dataset in RAM.

Shuffle semantics match shuffle_tfrecords_beam.py: every record is keyed by the
SHA1 of its bytes, records are hash-partitioned by that key into output shards,
and each shard is sorted by the key. The result is a deterministic, effectively
random global permutation.

Outputs (same layout as shuffle_tfrecords_beam.py):
  <output_pattern_prefix>-00000-of-000NN.tfrecord.gz ...
  <output_dataset_config_pbtxt>                     (optional)
  <dirname(output_pattern_prefix)>/example_info.json (copied if found)

Python deps on driver and executors: google-crc32c, boto3 (only for s3://).

Run locally:
  spark-submit --master 'local[*]' tools/shuffle_tfrecords_spark.py \
    --input_pattern_list="/tmp/training_set.with_label.tfrecord-?????-of-00016.gz" \
    --output_pattern_prefix="/tmp/out/training_set.with_label.shuffled" \
    --output_dataset_name="HG001" \
    --output_dataset_config_pbtxt="/tmp/out/training_set.dataset_config.pbtxt"

Run on EMR (see tools/aws_shuffle/README.md for cluster setup):
  spark-submit --deploy-mode cluster s3://BUCKET/shuffle_tfrecords_spark.py \
    --input_pattern_list="s3://BUCKET/training_set.with_label.tfrecord-?????-of-?????.gz" \
    --output_pattern_prefix="s3://BUCKET/shuffled/training_set.with_label.shuffled" \
    --output_dataset_name="HG001" \
    --output_dataset_config_pbtxt="s3://BUCKET/shuffled/training_set.dataset_config.pbtxt"
"""
# pylint: enable=line-too-long

import argparse
import collections
import fnmatch
import gzip
import hashlib
import io
import logging
import os
import re
import shutil
import struct
import textwrap

COMMENT_HEADER = """#
# --input_pattern_list={}
# --output_pattern_prefix={}
#
"""

# Default target size of one output shard, measured on compressed input bytes.
_DEFAULT_BYTES_PER_SHARD = 256 * 1024 * 1024
# S3 multipart part size. Must be >= 5 MiB.
_S3_PART_SIZE = 64 * 1024 * 1024


def parse_cmdline(argv=None):
  parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
  parser.add_argument(
      '--input_pattern_list',
      required=True,
      help='Comma-separated list of TFRecord filename patterns (s3:// or local).')
  parser.add_argument(
      '--output_pattern_prefix',
      required=True,
      help='Filename prefix for the output TFRecords.')
  parser.add_argument(
      '--output_dataset_config_pbtxt',
      help='Optional. If set, write a DeepVariantDatasetConfig pbtxt here.')
  parser.add_argument(
      '--output_dataset_name',
      help='Dataset name written into the pbtxt. Required with the pbtxt flag.')
  parser.add_argument(
      '--num_output_shards',
      type=int,
      default=0,
      help='Number of output shards (= shuffle partitions). 0 means derive it '
      'from total input size and --target_shard_mb.')
  parser.add_argument(
      '--target_shard_mb',
      type=int,
      default=_DEFAULT_BYTES_PER_SHARD // (1024 * 1024),
      help='Approximate compressed MB per output shard when '
      '--num_output_shards=0.')
  parser.add_argument(
      '--no_verify_crc',
      action='store_true',
      help='Skip CRC verification of input records.')
  args = parser.parse_args(argv)
  if args.output_dataset_config_pbtxt and not args.output_dataset_name:
    parser.error('Need to set --output_dataset_name.')
  return args


# ----------------------------------------------------------------------------
# Filesystem helpers: s3:// via boto3, everything else via the local filesystem.
# ----------------------------------------------------------------------------


def _split_s3(path):
  m = re.match(r'^s3a?://([^/]+)/?(.*)$', path)
  if not m:
    raise ValueError('Not an s3 path: {}'.format(path))
  return m.group(1), m.group(2)


def _is_s3(path):
  return path.startswith('s3://') or path.startswith('s3a://')


def _s3_client():
  import boto3  # pylint: disable=g-import-not-at-top
  return boto3.client('s3')


def glob_with_sizes(pattern):
  """Returns a sorted list of (path, size_bytes) matching a glob pattern."""
  if not _is_s3(pattern):
    import glob  # pylint: disable=g-import-not-at-top
    return sorted((p, os.path.getsize(p)) for p in glob.glob(pattern))
  bucket, key_pattern = _split_s3(pattern)
  # List from the longest literal prefix, then filter with fnmatch.
  literal_prefix = re.split(r'[*?\[]', key_pattern, maxsplit=1)[0]
  paginator = _s3_client().get_paginator('list_objects_v2')
  out = []
  for page in paginator.paginate(Bucket=bucket, Prefix=literal_prefix):
    for obj in page.get('Contents', []):
      if fnmatch.fnmatchcase(obj['Key'], key_pattern):
        out.append(('s3://{}/{}'.format(bucket, obj['Key']), obj['Size']))
  return sorted(out)


def open_read(path):
  if _is_s3(path):
    bucket, key = _split_s3(path)
    return _s3_client().get_object(Bucket=bucket, Key=key)['Body']
  return open(path, 'rb')


def copy_file(src, dst):
  if _is_s3(src) and _is_s3(dst):
    sb, sk = _split_s3(src)
    db, dk = _split_s3(dst)
    _s3_client().copy({'Bucket': sb, 'Key': sk}, db, dk)
  elif not _is_s3(src) and not _is_s3(dst):
    os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)
    shutil.copyfile(src, dst)
  else:
    with open_read(src) as f:
      write_bytes(dst, f.read())


def write_bytes(path, data):
  if _is_s3(path):
    bucket, key = _split_s3(path)
    _s3_client().put_object(Bucket=bucket, Key=key, Body=data)
  else:
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'wb') as f:
      f.write(data)


def exists(path):
  if _is_s3(path):
    bucket, key = _split_s3(path)
    resp = _s3_client().list_objects_v2(Bucket=bucket, Prefix=key, MaxKeys=1)
    return any(o['Key'] == key for o in resp.get('Contents', []))
  return os.path.exists(path)


class _S3MultipartWriter(io.RawIOBase):
  """Write-only file object that streams to S3 with multipart upload."""

  def __init__(self, path):
    super().__init__()
    self._bucket, self._key = _split_s3(path)
    self._client = _s3_client()
    self._buf = bytearray()
    self._parts = []
    self._upload_id = self._client.create_multipart_upload(
        Bucket=self._bucket, Key=self._key)['UploadId']

  def writable(self):
    return True

  def write(self, b):
    self._buf.extend(b)
    if len(self._buf) >= _S3_PART_SIZE:
      self._flush_part()
    return len(b)

  def _flush_part(self):
    n = len(self._parts) + 1
    resp = self._client.upload_part(
        Bucket=self._bucket, Key=self._key, UploadId=self._upload_id,
        PartNumber=n, Body=bytes(self._buf))
    self._parts.append({'PartNumber': n, 'ETag': resp['ETag']})
    self._buf = bytearray()

  def close(self):
    if self.closed:
      return
    try:
      if self._buf or not self._parts:
        self._flush_part()
      self._client.complete_multipart_upload(
          Bucket=self._bucket, Key=self._key, UploadId=self._upload_id,
          MultipartUpload={'Parts': self._parts})
    except Exception:
      self.abort()
      raise
    super().close()

  def abort(self):
    self._client.abort_multipart_upload(
        Bucket=self._bucket, Key=self._key, UploadId=self._upload_id)
    super().close()


def open_write(path):
  if _is_s3(path):
    return _S3MultipartWriter(path)
  os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
  return open(path, 'wb')


# ----------------------------------------------------------------------------
# TFRecord framing. Each record is:
#   uint64 length | uint32 masked_crc32c(length) | data | uint32 masked_crc32c(data)
# ----------------------------------------------------------------------------


def _masked_crc(data):
  import google_crc32c  # pylint: disable=g-import-not-at-top
  crc = google_crc32c.value(data)
  return (((crc >> 15) | (crc << 17)) + 0xa282ead8) & 0xffffffff


def read_tfrecords(path, verify_crc=True):
  """Yields raw record bytes from one (optionally gzipped) TFRecord file."""
  raw = open_read(path)
  f = gzip.GzipFile(fileobj=raw, mode='rb') if path.endswith('.gz') else raw
  try:
    while True:
      header = f.read(12)
      if not header:
        return
      if len(header) != 12:
        raise IOError('Truncated record header in {}'.format(path))
      length, length_crc = struct.unpack('<QI', header)
      if verify_crc and _masked_crc(header[:8]) != length_crc:
        raise IOError('Corrupt record length in {}'.format(path))
      data = f.read(length)
      footer = f.read(4)
      if len(data) != length or len(footer) != 4:
        raise IOError('Truncated record in {}'.format(path))
      if verify_crc and _masked_crc(data) != struct.unpack('<I', footer)[0]:
        raise IOError('Corrupt record data in {}'.format(path))
      yield data
  finally:
    f.close()
    raw.close()


def write_tfrecord(f, data):
  length = struct.pack('<Q', len(data))
  f.write(length)
  f.write(struct.pack('<I', _masked_crc(length)))
  f.write(data)
  f.write(struct.pack('<I', _masked_crc(data)))


# ----------------------------------------------------------------------------
# Minimal tf.Example decoder to read features['label'].int64_list.value[0]
# without importing TensorFlow on executors.
# ----------------------------------------------------------------------------


def _read_varint(buf, pos):
  result = shift = 0
  while True:
    b = buf[pos]
    pos += 1
    result |= (b & 0x7f) << shift
    if not b & 0x80:
      return result, pos
    shift += 7


def _iter_fields(buf):
  """Yields (field_number, wire_type, value) for a protobuf message."""
  pos, end = 0, len(buf)
  while pos < end:
    tag, pos = _read_varint(buf, pos)
    field, wire = tag >> 3, tag & 7
    if wire == 0:
      value, pos = _read_varint(buf, pos)
    elif wire == 1:
      value, pos = buf[pos:pos + 8], pos + 8
    elif wire == 2:
      n, pos = _read_varint(buf, pos)
      value, pos = buf[pos:pos + n], pos + n
    elif wire == 5:
      value, pos = buf[pos:pos + 4], pos + 4
    else:
      raise ValueError('Unsupported protobuf wire type {}'.format(wire))
    yield field, wire, value


def example_label(record):
  """Returns features['label'].int64_list.value[0], or None if absent."""
  buf = memoryview(record)
  for f1, _, features in _iter_fields(buf):  # Example.features = 1
    if f1 != 1:
      continue
    for f2, _, entry in _iter_fields(features):  # Features.feature = 1 (map)
      if f2 != 1:
        continue
      key = feature = None
      for f3, _, v in _iter_fields(entry):
        if f3 == 1:
          key = bytes(v)
        elif f3 == 2:
          feature = v
      if key != b'label' or feature is None:
        continue
      for f4, _, int64_list in _iter_fields(feature):  # Feature.int64_list = 3
        if f4 != 3:
          continue
        for f5, wire, v in _iter_fields(int64_list):  # Int64List.value = 1
          if f5 != 1:
            continue
          if wire == 2:  # packed
            v, _ = _read_varint(v, 0)
          return v - (1 << 64) if v >= (1 << 63) else v
  return None


# ----------------------------------------------------------------------------
# Pipeline.
# ----------------------------------------------------------------------------


def _shard_path(prefix, index, num_shards):
  return '{}-{:05d}-of-{:05d}.tfrecord.gz'.format(prefix, index, num_shards)


def _partition_of(key):
  return int.from_bytes(key[:8], 'big')


def make_config_string(name, tfrecord_path, num_examples):
  return textwrap.dedent("""
  name: "{}"
  tfrecord_path: "{}-?????-of-?????.tfrecord.gz"
  num_examples: {}
  """.format(name, tfrecord_path, num_examples))


def copy_example_info_json(input_pattern, output_pattern_prefix):
  matches = glob_with_sizes('{}*example_info.json'.format(input_pattern))
  if not matches:
    return
  dst = os.path.join(
      os.path.dirname(output_pattern_prefix), 'example_info.json')
  if not exists(dst):
    logging.info('Copying %s to %s', matches[0][0], dst)
    copy_file(matches[0][0], dst)


def main(argv=None):
  args = parse_cmdline(argv)
  from pyspark.sql import SparkSession  # pylint: disable=g-import-not-at-top

  spark = SparkSession.builder.appName('shuffle-tfrecords').getOrCreate()
  sc = spark.sparkContext

  patterns = [p.strip() for p in args.input_pattern_list.split(',') if p.strip()]
  copy_example_info_json(patterns[0], args.output_pattern_prefix)

  inputs = []
  for pattern in patterns:
    matched = glob_with_sizes(pattern)
    if not matched:
      raise ValueError('No files match {}'.format(pattern))
    inputs.extend(matched)
  total_bytes = sum(size for _, size in inputs)
  num_shards = args.num_output_shards or max(
      1, -(-total_bytes // (args.target_shard_mb * 1024 * 1024)))
  logging.info('Shuffling %d input files (%.2f GB compressed) into %d shards.',
               len(inputs), total_bytes / 1e9, num_shards)

  verify_crc = not args.no_verify_crc
  output_prefix = args.output_pattern_prefix

  def write_shard(index, sorted_pairs):
    labels = collections.Counter()
    writer = open_write(_shard_path(output_prefix, index, num_shards))
    try:
      with gzip.GzipFile(fileobj=writer, mode='wb') as gz:
        for _, record in sorted_pairs:
          write_tfrecord(gz, record)
          labels[example_label(record)] += 1
    except BaseException:
      if hasattr(writer, 'abort'):
        writer.abort()
      else:
        writer.close()
      raise
    writer.close()
    yield labels

  paths = [p for p, _ in inputs]
  label_counts = (
      sc.parallelize(paths, len(paths))
      .flatMap(lambda p: read_tfrecords(p, verify_crc))
      .map(lambda r: (hashlib.sha1(r).digest(), r))
      .repartitionAndSortWithinPartitions(
          numPartitions=num_shards, partitionFunc=_partition_of)
      .mapPartitionsWithIndex(write_shard)
      .reduce(lambda a, b: a + b))
  num_examples = sum(label_counts.values())
  logging.info('Wrote %d examples. Per label: %s', num_examples,
               dict(label_counts))

  if args.output_dataset_config_pbtxt:
    per_label = ''.join(
        '# class{}: {}\n'.format(k, v)
        for k, v in sorted(label_counts.items(), key=lambda kv: str(kv[0])))
    summary = ('# Generated by shuffle_tfrecords_spark.py\n' +
               COMMENT_HEADER.format(args.input_pattern_list, output_prefix) +
               per_label +
               make_config_string(args.output_dataset_name, output_prefix,
                                  num_examples))
    write_bytes(args.output_dataset_config_pbtxt, summary.encode('utf-8'))
  spark.stop()


if __name__ == '__main__':
  logging.basicConfig(level=logging.INFO)
  main()
