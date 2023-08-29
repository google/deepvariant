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
 */

#include "third_party/nucleus/core/statusor.h"

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

// Status tests.
TEST(StatusTest, OKStatusMatchesIsOK) {
  EXPECT_THAT(nucleus::Status(), IsOK());
}

TEST(StatusTest, FailedStatusMatchesIsNotOK) {
  EXPECT_THAT(Unknown("fail"), IsNotOK());
}

TEST(StatusTest, FailedStatusMatchesIsNotOKCode) {
  EXPECT_THAT(Unknown("fail"), IsNotOKWithCode(absl::StatusCode::kUnknown));
}

TEST(StatusTest, FailedStatusMatchesIsNotOKWithMessage) {
  EXPECT_THAT(Unknown("fail"), IsNotOKWithMessage("fail"));
}

TEST(StatusTest, FailedStatusMatchesIsNotOKWithCodeAndMessage) {
  EXPECT_THAT(Unknown("fail"),
              IsNotOKWithCodeAndMessage(absl::StatusCode::kUnknown, "fail"));
}

// StatusOr tests.
TEST(StatusOrTest, OKStatusMatchesIsOK) {
  StatusOr<int> status_or = StatusOr<int>(0);
  EXPECT_THAT(status_or, IsOK());
}

TEST(StatusOrTest, FailedStatusMatchesIsNotOK) {
  StatusOr<int> status_or = Unknown("fail");
  EXPECT_THAT(status_or, IsNotOK());
}

TEST(StatusOrTest, FailedStatusMatchesIsNotOKCode) {
  StatusOr<int> status_or = Unknown("fail");
  EXPECT_THAT(status_or, IsNotOKWithCode(absl::StatusCode::kUnknown));
}

TEST(StatusOrTest, FailedStatusMatchesIsNotOKWithMessage) {
  StatusOr<int> status_or = Unknown("fail");
  EXPECT_THAT(status_or, IsNotOKWithMessage("fail"));
}

TEST(StatusOrTest, FailedStatusMatchesIsNotOKWithCodeAndMessage) {
  StatusOr<int> status_or = Unknown("fail");
  EXPECT_THAT(status_or,
              IsNotOKWithCodeAndMessage(absl::StatusCode::kUnknown, "fail"));
}

}  // namespace nucleus
