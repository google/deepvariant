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
"""Common constants shared across DeepVariant's codebase.

This file is for very general constants in the code that end up needing to be
accessed in a variety of places, often in live code as well as throughout the
code in tests.
"""

from deepvariant.protos import deepvariant_pb2

# Default width [in basepairs] for our DeepVariant data tensor.
PILEUP_DEFAULT_WIDTH = 221

# Default height [in rows] for our DeepVariant data tensor.
PILEUP_DEFAULT_HEIGHT = 100

# Default channels if none are provided
PILEUP_DEFAULT_CHANNELS = [
    'read_base',
    'base_quality',
    'mapping_quality',
    'strand',
    'read_supports_variant',
    'base_differs_from_ref',
]

# Channels that include insert size.
PILEUP_CHANNELS_WITH_INSERT_SIZE = [
    'read_base',
    'base_quality',
    'mapping_quality',
    'strand',
    'read_supports_variant',
    'base_differs_from_ref',
    'insert_size',
]

# Not a default because it's hard-coded into the code.
PILEUP_NUM_CHANNELS = len(PILEUP_DEFAULT_CHANNELS)

# The dimensions of a pileup image tensor as height x width x rank.
PILEUP_DEFAULT_DIMS = [
    PILEUP_DEFAULT_HEIGHT,
    PILEUP_DEFAULT_WIDTH,
    PILEUP_NUM_CHANNELS,
]

# Number of classes represented in the data set. The three classes are
# homozygous reference (0), heterozygous (1) and homozygous alternative (2).
NUM_CLASSES = 3
NUM_DENOVO_CLASSES = 2

# Default sample name if no sample name is found from the BAM header.
DEFAULT_SAMPLE_NAME = 'default'

# Define all available channels, configured using --channel_list
CHANNELS = [
    'read_base',
    'base_quality',
    'mapping_quality',
    'strand',
    'read_supports_variant',
    'base_differs_from_ref',
    'haplotype',
    'allele_frequency',
    'diff_channels_alternate_allele_1',
    'diff_channels_alternate_allele_2',
    'read_mapping_percent',
    'avg_base_quality',
    'identity',
    'gap_compressed_identity',
    'gc_content',
    'is_homopolymer',
    'homopolymer_weighted',
    'blank',
    'insert_size',
    'base_channels_alternate_allele_1',
    'base_channels_alternate_allele_2',
    'mean_coverage',
]

ALT_ALIGNED_PILEUP_CHANNELS = [
    'base_channels_alternate_allele_1',
    'base_channels_alternate_allele_2',
    'diff_channels_alternate_allele_1',
    'diff_channels_alternate_allele_2',
]

# Create list of channels that can be used with --channel_list by removing
# channels specified using --alt_aligned_pileup.
USER_SET_CHANNELS = [
    x for x in CHANNELS if x not in ALT_ALIGNED_PILEUP_CHANNELS
]

CHANNEL_ENUM_TO_STRING = {
    deepvariant_pb2.CH_READ_BASE: 'read_base',
    deepvariant_pb2.CH_BASE_QUALITY: 'base_quality',
    deepvariant_pb2.CH_MAPPING_QUALITY: 'mapping_quality',
    deepvariant_pb2.CH_STRAND: 'strand',
    deepvariant_pb2.CH_READ_SUPPORTS_VARIANT: 'read_supports_variant',
    deepvariant_pb2.CH_BASE_DIFFERS_FROM_REF: 'base_differs_from_ref',
    deepvariant_pb2.CH_HAPLOTYPE_TAG: 'haplotype',
    deepvariant_pb2.CH_ALLELE_FREQUENCY: 'allele_frequency',
    deepvariant_pb2.CH_DIFF_CHANNELS_ALTERNATE_ALLELE_1: (
        'diff_channels_alternate_allele_1'
    ),
    deepvariant_pb2.CH_DIFF_CHANNELS_ALTERNATE_ALLELE_2: (
        'diff_channels_alternate_allele_2'
    ),
    deepvariant_pb2.CH_READ_MAPPING_PERCENT: 'read_mapping_percent',
    deepvariant_pb2.CH_AVG_BASE_QUALITY: 'avg_base_quality',
    deepvariant_pb2.CH_IDENTITY: 'identity',
    deepvariant_pb2.CH_GAP_COMPRESSED_IDENTITY: 'gap_compressed_identity',
    deepvariant_pb2.CH_GC_CONTENT: 'gc_content',
    deepvariant_pb2.CH_IS_HOMOPOLYMER: 'is_homopolymer',
    deepvariant_pb2.CH_HOMOPOLYMER_WEIGHTED: 'homopolymer_weighted',
    deepvariant_pb2.CH_BLANK: 'blank',
    deepvariant_pb2.CH_INSERT_SIZE: 'insert_size',
    deepvariant_pb2.CH_BASE_CHANNELS_ALTERNATE_ALLELE_1: (
        'base_channels_alternate_allele_1'
    ),
    deepvariant_pb2.CH_BASE_CHANNELS_ALTERNATE_ALLELE_2: (
        'base_channels_alternate_allele_2'
    ),
    deepvariant_pb2.CH_MEAN_COVERAGE: 'mean_coverage',
}

# Used only when phasing is on (phase_reads=true). It allows to set the
# region padding as a percantage over the region length. candidates are
# calculated over an extended region. Output examples are not affected by
# this value.
PHASE_READS_REGION_PADDING_PCT = 20
