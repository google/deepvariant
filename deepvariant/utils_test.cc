/*
 * Copyright 2017 Google LLC.
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
 */

#include "deepvariant/utils.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::EqualsProto;

TEST(UtilsTest, TestMakeAllele) {
  EXPECT_THAT(MakeAllele("A", AlleleType::REFERENCE, 1),
              EqualsProto("bases: \"A\" type: REFERENCE count: 1"));
  EXPECT_THAT(MakeAllele("AC", AlleleType::INSERTION, 10),
              EqualsProto("bases: \"AC\" type: INSERTION count: 10"));
  EXPECT_THAT(
      MakeAllele("AC", AlleleType::INSERTION, 10, true, 90, 5),
      EqualsProto("bases: \"AC\" type: INSERTION count: 10 is_low_quality: "
                  "true mapping_quality: 90 avg_base_quality: 5"));
}

TEST(UtilsTest, TestSimplifyRefAlt) {
  EXPECT_EQ(SimplifyRefAlt("CAA", "CA"), "CA->C");
  EXPECT_EQ(SimplifyRefAlt("CA", "C"), "CA->C");
  EXPECT_EQ(SimplifyRefAlt("ATGTG", "ATGTGTGTGTGTG"), "A->ATGTGTGTG");
  // This input will likely not happen the way we currently use this function.
  // But I want to make sure the output is not empty.
  EXPECT_EQ(SimplifyRefAlt("C", "C"), "C->C");
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
