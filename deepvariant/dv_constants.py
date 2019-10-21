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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# Default width [in basepairs] for our DeepVariant data tensor.
PILEUP_DEFAULT_WIDTH = 221

# Default height [in rows] for our DeepVariant data tensor.
PILEUP_DEFAULT_HEIGHT = 100

# Not a default because it's hard-coded into the code.
PILEUP_NUM_CHANNELS = 6

# The dimensions of a pileup image tensor as height x width x rank.
PILEUP_DEFAULT_DIMS = [
    PILEUP_DEFAULT_HEIGHT, PILEUP_DEFAULT_WIDTH, PILEUP_NUM_CHANNELS
]

# Number of classes represented in the data set. The three classes are
# homozygous reference (0), heterozygous (1) and homozygous alternative (2).
NUM_CLASSES = 3

# Default sample name if no sample name is found from the BAM header.
DEFAULT_SAMPLE_NAME = 'default'
