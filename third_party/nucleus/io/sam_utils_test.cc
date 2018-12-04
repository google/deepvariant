/*
 * Copyright 2018 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "third_party/nucleus/io/sam_utils.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "htslib/sam.h"
#include "third_party/nucleus/protos/cigar.pb.h"

namespace nucleus {

using genomics::v1::CigarUnit;
using genomics::v1::CigarUnit_Operation_Operation_MAX;
using genomics::v1::CigarUnit_Operation_Operation_MIN;

TEST(SamUtilsTest, Conversion) {
  for (int i = CigarUnit_Operation_Operation_MIN;
       i <= CigarUnit_Operation_Operation_MAX; ++i) {
    CigarUnit::Operation op = static_cast<CigarUnit::Operation>(i);
    EXPECT_EQ(op, kHtslibCigarToProto[kProtoToHtslibCigar[i]]);
  }
}

TEST(SamUtilsTest, ProtoConversion) {
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::ALIGNMENT_MATCH], BAM_CMATCH);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::INSERT], BAM_CINS);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::DELETE], BAM_CDEL);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::SKIP], BAM_CREF_SKIP);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::CLIP_SOFT], BAM_CSOFT_CLIP);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::CLIP_HARD], BAM_CHARD_CLIP);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::PAD], BAM_CPAD);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::SEQUENCE_MATCH], BAM_CEQUAL);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::SEQUENCE_MISMATCH], BAM_CDIFF);
  EXPECT_EQ(kProtoToHtslibCigar[CigarUnit::OPERATION_UNSPECIFIED], BAM_CBACK);
}

}  // namespace nucleus
