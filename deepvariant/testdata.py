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

DEEPVARIANT_DATADIR = ''


def deepvariant_testdata(filename):
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
      os.path.join('deepvariant/testdata', filename), DEEPVARIANT_DATADIR
  )


CHR20_FASTA = None
CHR20_BAM = None
CHR20_BAM_FIRST_HALF = None
CHR20_BAM_SECOND_HALF = None
NOCHR20_BAM = None
CHR20_CRAM = None
GOLDEN_TRAINING_EXAMPLES = None
GOLDEN_CALLING_CANDIDATES = None
GOLDEN_CANDIDATE_POSITIONS = None
GOLDEN_CALLING_EXAMPLES = None
CONFIDENT_REGIONS_BED = None
TRUTH_VARIANTS_VCF = None
TRUTH_VARIANTS_VCF_WITH_TYPES = None
GOLDEN_POSTPROCESS_INPUT = None
GOLDEN_POSTPROCESS_INPUT_SHARDED = None
GOLDEN_POSTPROCESS_OUTPUT = None
GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY = None
GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED = None
GOLDEN_POSTPROCESS_GVCF_INPUT = None
GOLDEN_POSTPROCESS_GVCF_OUTPUT = None
GOLDEN_POSTPROCESS_GVCF_OUTPUT_COMPRESSED = None
GOLDEN_MAKE_EXAMPLES_RUN_INFO = None
WS_ALLELE_COUNT_LINEAR_MODEL = None
WS_ALLELE_COUNT_LINEAR_MODEL_PCKL = None
WS_VARIANT_READS_THRESHOLD_MODEL = None
GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT = None
GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT = None

N_GOLDEN_TRAINING_EXAMPLES = 49
N_GOLDEN_CALLING_EXAMPLES = 84

# For CustomizedClassesVariantLabeler:
CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES = None

# For VcfCandidateImporter:
GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES = None
GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES = None
VCF_CANDIDATE_IMPORTER_VARIANTS = None

# For alt-aligned pileups:
ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES = None
ALT_ALIGNED_ROWS_EXAMPLES = None
RUNTIME_BY_REGION = None
RUNTIME_BY_REGION_SHARDED = None

# For allele frequency:
VCF_WITH_ALLELE_FREQUENCIES = None
GRCH38_FASTA = None
AF_VCF_CHR20 = None
AF_VCF_CHR21 = None
AF_VCF_CHR20_21_WILDCARD = None
AF_VCF_CHR20_AND_21 = None
GRCH38_CHR20_AND_21_BAM = None
GOLDEN_ALLELE_FREQUENCY_EXAMPLES = None



def init():
  """Initialize global variables from flag values."""
  global CHR20_FASTA
  global CHR20_BAM
  global CHR20_BAM_FIRST_HALF
  global CHR20_BAM_SECOND_HALF
  global NOCHR20_BAM
  global CHR20_CRAM
  global GOLDEN_TRAINING_EXAMPLES
  global GOLDEN_CALLING_CANDIDATES
  global GOLDEN_CANDIDATE_POSITIONS
  global GOLDEN_CALLING_EXAMPLES
  global CONFIDENT_REGIONS_BED
  global TRUTH_VARIANTS_VCF
  global TRUTH_VARIANTS_VCF_WITH_TYPES
  global GOLDEN_POSTPROCESS_INPUT
  global GOLDEN_POSTPROCESS_INPUT_SHARDED
  global GOLDEN_POSTPROCESS_OUTPUT
  global GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY
  global GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED
  global GOLDEN_POSTPROCESS_GVCF_INPUT
  global GOLDEN_POSTPROCESS_GVCF_OUTPUT
  global GOLDEN_POSTPROCESS_GVCF_OUTPUT_COMPRESSED
  global GOLDEN_MAKE_EXAMPLES_RUN_INFO
  global WS_ALLELE_COUNT_LINEAR_MODEL
  global WS_ALLELE_COUNT_LINEAR_MODEL_PCKL
  global WS_VARIANT_READS_THRESHOLD_MODEL
  global GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT
  global GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT

  CHR20_FASTA = deepvariant_testdata('input/ucsc.hg19.chr20.unittest.fasta.gz')
  CHR20_BAM = deepvariant_testdata('input/NA12878_S1.chr20.10_10p1mb.bam')
  # # Here is how "NA12878_S1.chr20.10_10p1mb.first_half.bam"
  # # and "NA12878_S1.chr20.10_10p1mb.second_half.bam" are split
  # # from NA12878_S1.chr20.10_10p1mb.bam.
  # READS_FIRST_HALF=${TESTDATA_DIR}/NA12878_S1.chr20.10_10p1mb.first_half.bam
  # READS_SECOND_HALF=${TESTDATA_DIR}/NA12878_S1.chr20.10_10p1mb.second_half.bam
  # READS=${TESTDATA_DIR}/NA12878_S1.chr20.10_10p1mb.bam
  # samtools view -H ${READS} > /tmp/f1.sam
  # cp /tmp/f1.sam /tmp/f2.sam
  # # Because ${READS} has total of 52035 lines, we split in roughly half.
  # samtools view ${READS} | head -26000 >> /tmp/f1.sam
  # samtools view ${READS} | tail -26035 >> /tmp/f2.sam
  # samtools view -S -b /tmp/f1.sam > ${READS_FIRST_HALF}
  # samtools view -S -b /tmp/f2.sam > ${READS_SECOND_HALF}
  # samtools index ${READS_FIRST_HALF}
  # samtools index ${READS_SECOND_HALF}
  CHR20_BAM_FIRST_HALF = deepvariant_testdata(
      'input/NA12878_S1.chr20.10_10p1mb.first_half.bam'
  )
  CHR20_BAM_SECOND_HALF = deepvariant_testdata(
      'input/NA12878_S1.chr20.10_10p1mb.second_half.bam'
  )
  # # Here is how the "HG002_NIST_150bp_downsampled_30x.chr20.10_10p1mb.bam"
  # # file was created.
  # samtools view -hb HG002_NIST_150bp_downsampled_30x.bam \
  #     20:10,000,000-10,100,000 \
  #     > HG002_NIST_150bp_downsampled_30x.chr20.10_10p1mb.bam
  # samtools index HG002_NIST_150bp_downsampled_30x.chr20.10_10p1mb.bam
  NOCHR20_BAM = deepvariant_testdata(
      'input/HG002_NIST_150bp_downsampled_30x.chr20.10_10p1mb.bam'
  )
  CHR20_CRAM = deepvariant_testdata('input/NA12878_S1.chr20.10_10p1mb.cram')
  GOLDEN_TRAINING_EXAMPLES = deepvariant_testdata(
      'golden.training_examples.tfrecord.gz'
  )
  GOLDEN_CALLING_CANDIDATES = deepvariant_testdata(
      'golden.calling_examples.tfrecord.gz'
  )
  GOLDEN_CANDIDATE_POSITIONS = deepvariant_testdata(
      'golden.candidate_positions'
  )
  GOLDEN_CALLING_EXAMPLES = deepvariant_testdata(
      'golden.calling_examples.tfrecord.gz'
  )
  CONFIDENT_REGIONS_BED = deepvariant_testdata(
      'input/test_nist.b37_chr20_100kbp_at_10mb.bed'
  )
  TRUTH_VARIANTS_VCF = deepvariant_testdata(
      'input/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz'
  )
  TRUTH_VARIANTS_VCF_WITH_TYPES = deepvariant_testdata(
      'input/with_types.test_nist.b37_chr20_4kbp_at_10mb.vcf.gz'
  )
  GOLDEN_POSTPROCESS_INPUT = deepvariant_testdata(
      'golden.postprocess_single_site_input.tfrecord.gz'
  )
  GOLDEN_POSTPROCESS_INPUT_SHARDED = deepvariant_testdata(
      'golden.postprocess_single_site_input-00000-of-00001.tfrecord.gz'
  )
  GOLDEN_POSTPROCESS_OUTPUT = deepvariant_testdata(
      'golden.postprocess_single_site_output.vcf'
  )
  GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY = deepvariant_testdata(
      'golden.postprocess_single_site_output.pass_only.vcf'
  )
  GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED = deepvariant_testdata(
      'golden.postprocess_single_site_output.vcf.gz'
  )
  GOLDEN_POSTPROCESS_GVCF_INPUT = deepvariant_testdata(
      'golden.postprocess_gvcf_input.tfrecord.gz'
  )
  GOLDEN_POSTPROCESS_GVCF_OUTPUT = deepvariant_testdata(
      'golden.postprocess_gvcf_output.g.vcf'
  )
  GOLDEN_MAKE_EXAMPLES_RUN_INFO = deepvariant_testdata(
      'golden.training_examples.tfrecord.gz.run_info.pbtxt'
  )
  WS_ALLELE_COUNT_LINEAR_MODEL = deepvariant_testdata(
      'obsolete/window_selector_allele_count_linear.pbtxt'
  )
  WS_ALLELE_COUNT_LINEAR_MODEL_PCKL = deepvariant_testdata(
      'obsolete/window_selector_allele_count_linear.pckl'
  )
  WS_VARIANT_READS_THRESHOLD_MODEL = deepvariant_testdata(
      'obsolete/window_selector_variant_read_threshold.pbtxt'
  )
  GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT = deepvariant_testdata(
      'golden.vcf_candidate_importer_postprocess_single_site_input.tfrecord.gz'
  )
  GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT = deepvariant_testdata(
      'golden.vcf_candidate_importer_postprocess_single_site_output.vcf'
  )

  # For CustomizedClassesVariantLabeler:
  global CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES
  CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES = deepvariant_testdata(
      'customized_classes.golden.training_examples.tfrecord.gz'
  )

  # For VcfCandidateImporter:
  global GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES
  global GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES
  global VCF_CANDIDATE_IMPORTER_VARIANTS
  GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES = deepvariant_testdata(
      'golden.vcf_candidate_importer.training_examples.tfrecord.gz'
  )
  GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES = deepvariant_testdata(
      'golden.vcf_candidate_importer_calling_examples.tfrecord'
  )
  VCF_CANDIDATE_IMPORTER_VARIANTS = deepvariant_testdata(
      'input/vcf_candidate_importer.indels.chr20.vcf.gz'
  )

  # For alt-aligned pileups:
  global ALT_ALIGNED_ROWS_EXAMPLES
  global ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES
  ALT_ALIGNED_ROWS_EXAMPLES = deepvariant_testdata(
      'golden.alt_aligned_pileup_rows_examples.tfrecord.gz'
  )
  ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES = deepvariant_testdata(
      'golden.alt_aligned_pileup_diff_channels_examples.tfrecord.gz'
  )

  # For runtime-by-region in make_examples:
  global RUNTIME_BY_REGION
  global RUNTIME_BY_REGION_SHARDED
  RUNTIME_BY_REGION = deepvariant_testdata('input/make_examples_runtime.tsv')
  RUNTIME_BY_REGION_SHARDED = deepvariant_testdata(
      'input/make_examples_runtime@2.tsv'
  )

  # For allele_frequency with GRCh38:
  global VCF_WITH_ALLELE_FREQUENCIES
  global GRCH38_FASTA
  global AF_VCF_CHR20
  global AF_VCF_CHR21
  global AF_VCF_CHR20_21_WILDCARD
  global AF_VCF_CHR20_AND_21
  global GRCH38_CHR20_AND_21_BAM
  global GOLDEN_ALLELE_FREQUENCY_EXAMPLES
  VCF_WITH_ALLELE_FREQUENCIES = deepvariant_testdata(
      'input/allele_frequencies_vcf.vcf.gz'
  )

  # Fasta filtered to regions: chr20:1-10000000 and chr21:1-10000000.
  GRCH38_FASTA = deepvariant_testdata('input/grch38.chr20_and_21_10M.fa.gz')
  # VCFs filtered to chr20:1-100000 and chr21:5100000-5200000.
  AF_VCF_CHR20 = deepvariant_testdata('input/cohort-chr20_100k.vcf.gz')
  AF_VCF_CHR21 = deepvariant_testdata('input/cohort-chr21_100k.vcf.gz')
  AF_VCF_CHR20_AND_21 = deepvariant_testdata(
      'input/cohort-chr20_and_chr21_100k.vcf.gz'
  )
  AF_VCF_CHR20_21_WILDCARD = deepvariant_testdata(
      'input/cohort-chr2?_100k.vcf.gz'
  )
  # This bam filtered to regions: chr20:61001-62000 and chr21:5114000-5114999
  # and header is edited with the following to match the GRCH38_FASTA:
  # @SQ     SN:chr20        LN:10000000
  # @SQ     SN:chr21        LN:10000000
  GRCH38_CHR20_AND_21_BAM = deepvariant_testdata(
      'input/grch38_1k_subset_chr20_and_chr21.bam'
  )
  GOLDEN_ALLELE_FREQUENCY_EXAMPLES = deepvariant_testdata(
      'golden.allele_frequency_examples.tfrecord.gz'
  )

