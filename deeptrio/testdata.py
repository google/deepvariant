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
"""Utilities to help with testing DeepVariant code."""

import os



from third_party.nucleus.testing import test_utils as nucleus_test_utils

GENOMICS_DIR = 'learning/genomics'


def deeptrio_testdata(filename):
  """Gets the path to filename in genomics/deepvariant/testdata.

  These paths are only known at runtime, after flag parsing
  has occurred.

  Args:
    filename: The name of a testdata file in the core genomics testdata
      directory. For example, if you have a test file in
      "learning/genomics/deepvariant/testdata/foo.txt", filename should be
      "foo.txt" to get a path to it.

  Returns:
    The absolute path to a testdata file.
  """
  return nucleus_test_utils.genomics_testdata(
      os.path.join('deeptrio/testdata', filename), GENOMICS_DIR
  )


CHR20_FASTA = None
HG001_CHR20_BAM = None
NA12891_CHR20_BAM = None
NA12892_CHR20_BAM = None
GOLDEN_TRAINING_EXAMPLES = None
GOLDEN_CALLING_CANDIDATES = None
GOLDEN_CANDIDATE_POSITIONS = None
GOLDEN_CALLING_EXAMPLES = None
CONFIDENT_REGIONS_BED = None
TRUTH_VARIANTS_VCF = None
TRUTH_VARIANTS_VCF_WITH_TYPES = None
GOLDEN_POSTPROCESS_INPUT = None
GOLDEN_POSTPROCESS_OUTPUT = None
GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED = None
GOLDEN_POSTPROCESS_GVCF_INPUT = None
GOLDEN_POSTPROCESS_GVCF_OUTPUT = None
GOLDEN_POSTPROCESS_GVCF_OUTPUT_COMPRESSED = None
GOLDEN_MAKE_EXAMPLES_RUN_INFO = None
WS_ALLELE_COUNT_LINEAR_MODEL = None
WS_ALLELE_COUNT_LINEAR_MODEL_PCKL = None
WS_VARIANT_READS_THRESHOLD_MODEL = None
# Test data for ONT
GRCH38_CHR0_FASTA = None
ONT_HG002_BAM = None
ONT_HG003_BAM = None
ONT_HG004_BAM = None
HG002_HIGH_CONFIDENCE_VCF = None
HG002_HIGH_CONFIDENCE_BED = None
HG002_DENOVO_BED = None
GOLDEN_ONT_MAKE_EXAMPLES_OUTPUT = None
GOLDEN_ONT_DENOVO_MAKE_EXAMPLES_OUTPUT = None

ONT_N_GOLDEN_TRAINING_EXAMPLES = 167
N_GOLDEN_TRAINING_EXAMPLES = 50
N_GOLDEN_CALLING_EXAMPLES = 103

CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES = None
ALT_ALIGNED_PILEUP_GOLDEN_TRAINING_EXAMPLES = None
GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES = None
GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES_CHILD = None


def init():
  """Initialize global variables from flag values."""
  global CHR20_FASTA
  global HG001_CHR20_BAM
  global NA12891_CHR20_BAM
  global NA12892_CHR20_BAM
  global GOLDEN_TRAINING_EXAMPLES
  global GOLDEN_CANDIDATE_POSITIONS
  global GOLDEN_CALLING_CANDIDATES
  global GOLDEN_CALLING_EXAMPLES
  global CONFIDENT_REGIONS_BED
  global TRUTH_VARIANTS_VCF
  global TRUTH_VARIANTS_VCF_WITH_TYPES
  global GOLDEN_POSTPROCESS_INPUT
  global GOLDEN_POSTPROCESS_OUTPUT
  global GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED
  global GOLDEN_POSTPROCESS_GVCF_INPUT
  global GOLDEN_POSTPROCESS_GVCF_OUTPUT
  global GOLDEN_POSTPROCESS_GVCF_OUTPUT_COMPRESSED
  global GOLDEN_MAKE_EXAMPLES_RUN_INFO
  global WS_ALLELE_COUNT_LINEAR_MODEL
  global WS_ALLELE_COUNT_LINEAR_MODEL_PCKL
  global WS_VARIANT_READS_THRESHOLD_MODEL
  global GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES
  global GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES_CHILD

  global GRCH38_CHR0_FASTA
  global ONT_HG002_BAM
  global ONT_HG003_BAM
  global ONT_HG004_BAM
  global HG002_HIGH_CONFIDENCE_VCF
  global HG002_HIGH_CONFIDENCE_BED
  global HG002_DENOVO_BED
  global GOLDEN_ONT_MAKE_EXAMPLES_OUTPUT
  global GOLDEN_ONT_DENOVO_MAKE_EXAMPLES_OUTPUT

  CHR20_FASTA = deeptrio_testdata('input/hs37d5.chr20.fa.gz')
  HG001_CHR20_BAM = deeptrio_testdata('input/HG001.chr20.10_10p1mb_sorted.bam')
  NA12891_CHR20_BAM = deeptrio_testdata(
      'input/NA12891.chr20.10_10p1mb_sorted.bam'
  )
  NA12892_CHR20_BAM = deeptrio_testdata(
      'input/NA12892.chr20.10_10p1mb_sorted.bam'
  )

  GOLDEN_TRAINING_EXAMPLES = deeptrio_testdata(
      'golden.training_examples.tfrecord.gz'
  )
  GOLDEN_CANDIDATE_POSITIONS = deeptrio_testdata(
      'golden_child.candidate_positions'
  )
  GOLDEN_CALLING_CANDIDATES = deeptrio_testdata(
      'golden_child.calling_examples.tfrecord.gz'
  )
  GOLDEN_CALLING_EXAMPLES = deeptrio_testdata(
      'golden_child.calling_examples.tfrecord.gz'
  )
  CONFIDENT_REGIONS_BED = deeptrio_testdata(
      'input/test_giab.b37_chr20_100kbp_at_10mb.bed'
  )
  TRUTH_VARIANTS_VCF = deeptrio_testdata(
      'input/HG001_chr20_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
  )
  TRUTH_VARIANTS_VCF_WITH_TYPES = deeptrio_testdata(
      'input/with_types.test_nist.b37_chr20_4kbp_at_10mb.vcf.gz'
  )
  GOLDEN_POSTPROCESS_INPUT = deeptrio_testdata(
      'golden.postprocess_single_site_input.tfrecord.gz'
  )
  GOLDEN_POSTPROCESS_OUTPUT = deeptrio_testdata(
      'golden.postprocess_single_site_output.vcf'
  )
  GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED = deeptrio_testdata(
      'golden.postprocess_single_site_output.vcf.gz'
  )
  GOLDEN_POSTPROCESS_GVCF_INPUT = deeptrio_testdata(
      'golden_child.postprocess_gvcf_input.tfrecord.gz'
  )
  GOLDEN_POSTPROCESS_GVCF_OUTPUT = deeptrio_testdata(
      'golden.postprocess_gvcf_output.g.vcf'
  )
  GOLDEN_MAKE_EXAMPLES_RUN_INFO = deeptrio_testdata(
      'golden.training_examples.tfrecord.gz.run_info.pbtxt'
  )
  WS_ALLELE_COUNT_LINEAR_MODEL = deeptrio_testdata(
      'window_selector_allele_count_linear.pbtxt'
  )
  WS_ALLELE_COUNT_LINEAR_MODEL_PCKL = deeptrio_testdata(
      'window_selector_allele_count_linear.pckl'
  )
  WS_VARIANT_READS_THRESHOLD_MODEL = deeptrio_testdata(
      'window_selector_variant_read_threshold.pbtxt'
  )

  # For oxford nanopore
  GRCH38_CHR0_FASTA = deeptrio_testdata(
      'input/grch38.chr20_5050000_5075000.masked.fa.gz'
  )
  ONT_HG002_BAM = deeptrio_testdata('input/HG002_R10_chr20_5050000_5075000.bam')
  ONT_HG003_BAM = deeptrio_testdata('input/HG003_R10_chr20_5050000_5075000.bam')
  ONT_HG004_BAM = deeptrio_testdata('input/HG004_R10_chr20_5050000_5075000.bam')
  HG002_HIGH_CONFIDENCE_VCF = deeptrio_testdata(
      'input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz'
  )
  HG002_HIGH_CONFIDENCE_BED = deeptrio_testdata(
      'input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.bed'
  )
  HG002_DENOVO_BED = deeptrio_testdata(
      'input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.denovo_regions.bed'
  )
  GOLDEN_ONT_MAKE_EXAMPLES_OUTPUT = deeptrio_testdata(
      'HG002_ONT_deeptrio.examples.tfrecord.gz'
  )
  GOLDEN_ONT_DENOVO_MAKE_EXAMPLES_OUTPUT = deeptrio_testdata(
      'HG002_ONT_deeptrio.denovo.examples.tfrecord.gz'
  )

  # For CustomizedClassesVariantLabeler.
  global CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES
  CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES = deeptrio_testdata(
      'customized_classes.golden.training_examples.tfrecord.gz'
  )

  # For alt-aligned pileups
  global ALT_ALIGNED_PILEUP_GOLDEN_TRAINING_EXAMPLES
  ALT_ALIGNED_PILEUP_GOLDEN_TRAINING_EXAMPLES = deeptrio_testdata(
      'alt_aligned_pileup.golden.training_examples.tfrecord.gz'
  )

  GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES = deeptrio_testdata(
      'golden.vcf_candidate_importer.training_examples.tfrecord.gz'
  )
  GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES_CHILD = deeptrio_testdata(
      'golden_child.vcf_candidate_importer.calling_examples.tfrecord.gz'
  )
