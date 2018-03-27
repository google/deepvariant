r"""Creates tf.Example protos for learning the error process of a sequencer.

This program reads in aligned NGS reads from --bam and emits TFRecord of
tf.Examples containing the observed bases and quality scores for each read as
well as the true genome sequence for that read, which are determined from the
reference genome (--ref) and the known variants for the sample sequenced in the
BAM (--vcf). For example, suppose we have a read in --bam:

  sequence: ACGT
  aligned at: chr20:10-15
  cigar: 4M

And that the reference genome at chr20:10-15 is ACCT. This program first checks
that --vcf file doesn't have any genetic variants in the region chr20:10-15. If
there's no variants there, the program emits a tf.Example containing:

    read_sequence - the bases of the read ('ACGT')
    true_sequence - the reference genome bases ('ACCT')
    as well as read_name, cigar, and read_qualities

When run over a BAM file, we generate one tf.Example for each read that
satisfies our read requirements (e.g., aligned, not a duplicate), that spans a
region of the genome without genetic variants in --vcf, and spans reference
bases containing only A, C, G, or T bases.

The emitted examples are suitable (but not ideal, see below) for training a
TensorFlow model to predict the true sequence for a read. This learned model can
be considered a model for the error process of the sequencer. It can be used to
to "error correct" a read, or to compute a distribution of possible true
sequences for a given read.

NOTE: we have not trained a TensorFlow model based on this program. It is an
interesting technical / scientific question what kind of model architecture
would perform well here, and how accurately it can learn to identify errors in
the input reads. There are many extensions to this example program that would be
important to add when attempting to training a production-grade error model,
though (see next for one example).

NOTE: This program doesn't filter out reads that aren't in the confidently
determined region of the genome. A more sophisticated extension would also
accept a confident regions BED file and use ranges.RangeSet to only include
reads fully within the confident regions. This is left as an exercise for the
user.

## Prerequisites

a) Grab a slice of the BAM file via samtools. We don't need to index it because
   we are only iterating over the records, and don't need the query()
   functionality.

samtools view -h \
  ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam \
  20:10,000,000-10,100,000 \
  -o /tmp/NA12878_sliced.bam

b) Get the Genome in a Bottle truth set for NA12878 and index it with tabix.

wget \
  ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
  -O /tmp/NA12878_calls.vcf.gz

tabix /tmp/NA12878_calls.vcf.gz

c) Get the reference genome FASTA file from the DeepVariant testdata GCP bucket
   and index it with samtools.

gsutil cp gs://deepvariant/case-study-testdata/hs37d5.fa.gz /tmp/hs37d5.fa.gz
samtools faidx /tmp/hs37d5.fa.gz

## Build and run

blaze build -c opt //third_party/nucleus/examples:ngs_errors

./blaze-bin/third_party/nucleus/examples/ngs_errors \
  --alsologtostderr \
  --ref /tmp/hs37d5.fa.gz \
  --vcf /tmp/NA12878_calls.vcf.gz \
  --bam /tmp/NA12878_sliced.bam \
  --examples_out /tmp/examples.tfrecord
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from absl import app
from tensorflow import flags
from absl import logging

from deepvariant.util.io import fasta
from deepvariant.util.io import genomics_writer
from deepvariant.util.io import sam
from deepvariant.util.io import vcf
from deepvariant.util.genomics import sam_pb2
from deepvariant.util import cigar
from deepvariant.util import ranges
from deepvariant.util import utils
from tensorflow.core.example import example_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'bam',
    None,
    'Path to an aligned SAM/BAM file (not necessarily indexed) containing '
    'reads to convert to tf.Examples.')
flags.DEFINE_string(
    'ref',
    None,
    'Path to an indexed FASTA file to use to get the reference sequence for '
    'the reads in --bam. Should be the same reference genome used to align the '
    'SAM/BAM file provided by --bam.')
flags.DEFINE_string(
    'vcf',
    None,
    'Path to an indexed VCF file containing variant calls for the sample in '
    '--bam. These variants are used to eliminate regions where the sample '
    'differs from the reference genome.')
flags.DEFINE_string(
    'examples_out',
    None,
    'Path where we will write our TFRecord of tf.Examples.')

_ALLOWED_CIGAR_OPS = frozenset([cigar.CHAR_TO_CIGAR_OPS[op] for op in 'MX=ID'])


def is_usable_training_example(read, variants, ref_bases):
  """Returns True if we can use read to make a training example.

  Args:
    read: nucleus.Read proto.
    variants: list[nucleus.Variant] protos. A list of variants overlapping read.
    ref_bases: str. The reference bases for read.

  Returns:
    True if read can be used to construct a high-quality training example, False
    otherwise.
  """
  if variants:
    # We cannot be certain that the true sequence of read is the bases from the
    # reference genome because there are overlapping variant calls, indicating
    # that the sample may differ from the reference sequence under this read.
    return False
  if not all(cigarunit.operation in _ALLOWED_CIGAR_OPS
             for cigarunit in read.alignment.cigar):
    # The read's alignment is too complex for us to use, so return False.
    return False
  if any(base not in 'ACGT' for base in ref_bases):
    # The reference bases aren't just one of the ACGT option, so return False.
    return False
  return True


def make_example(read, ref_bases):
  """Create a tf.Example for read and ref_bases.

  Args:
    read: nucleus.Read proto with cigar, fragment_name, and aligned_sequence.
    variants: list[nucleus.Variant] protos. A list of variants overlapping read.
    ref_bases: str. The reference bases for read.

  Returns:
    A tf.Example protobuf with the following features:
      read_name - for debugging convenience
      cigar - the cigar string of the read
      read_sequence - the bases observed by the instrument
      read_qualities - the quality scores emitted by the instrument, as
                       phred-scaled integer values.
      true_sequence - the "true" bases that should have been observed for this
                      read, as extracted from the reference genome.
  """
  # Create our example proto.
  example = example_pb2.Example()

  # Set the features in our Example.
  # Note that the str(...) calls are necessary because proto string fields are
  # unicode objects and we can only add bytes to the bytes_list.
  features = example.features
  features.feature['read_name'].bytes_list.value.append(str(read.fragment_name))
  features.feature['cigar'].bytes_list.value.append(
      cigar.format_cigar_units(read.alignment.cigar))
  features.feature['read_sequence'].bytes_list.value.append(
      str(read.aligned_sequence))
  features.feature['read_qualities'].int64_list.value.extend(
      read.aligned_quality)
  features.feature['true_sequence'].bytes_list.value.append(ref_bases)

  return example


def main(argv):
  del argv  # Unused.

  # Create a ref_reader backed by FLAGS.ref.
  ref_reader = fasta.RefFastaReader(FLAGS.ref)

  # Create a vcf_reader backed by FLAGS.vcf.
  vcf_reader = vcf.VcfReader(FLAGS.vcf)

  # Create a sam_reader backed by FLAGS.bam. Provide an empty ReadRequirements
  # proto to the reader so it enables standard filtering based on the default
  # values of ReadRequirements. Also explicitly allow the reader to access an
  # unindexed BAM, so only the iterate() function is enabled.
  read_requirements = sam_pb2.ReadRequirements()
  sam_reader = sam.SamReader(
      FLAGS.bam, use_index=False, read_requirements=read_requirements)

  # Create our TFRecordWriter where we'll send our tf.Examples.
  examples_out = genomics_writer.TFRecordWriter(FLAGS.examples_out)

  # All our readers and writers are context managers, so use the `with`
  # construct to open all of the inputs/outputs and close them when we are done
  # looping over our reads.
  n_examples = 0
  with ref_reader, vcf_reader, sam_reader, examples_out:
    # loop over the reads in our BAM file:
    for i, read in enumerate(sam_reader.iterate(), start=1):
      # Get the Range proto describing the chrom/start/stop spanned by our read.
      read_range = utils.read_range(read)

      # Get all of the variants that overlap our read range.
      variants = list(vcf_reader.query(read_range))

      # Get the reference bases spanned by our read.
      ref_bases = ref_reader.query(read_range)

      # Check that we can use our read for generating a example.
      if is_usable_training_example(read, variants, ref_bases):
        n_examples += 1

        # Convert read and ref_bases to a tf.Example with make_example.
        example = make_example(read, ref_bases)

        # And write it out to our TFRecord output file.
        examples_out.write(example)

        # Do a bit of convenient logging. This is very verbose if we convert a
        # lot of reads...
        logging.info(
            ('Added an example for read %s (span=%s) with cigar %s [%d added '
             'of %d total reads]'), read.fragment_name,
            ranges.to_literal(read_range),
            cigar.format_cigar_units(read.alignment.cigar), n_examples, i)


if __name__ == '__main__':
  flags.mark_flags_as_required(['ref', 'vcf', 'bam', 'examples_out'])
  app.run(main)
