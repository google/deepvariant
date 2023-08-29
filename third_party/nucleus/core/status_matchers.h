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

#ifndef THIRD_PARTY_NUCLEUS_VENDOR_STATUS_MATCHERS_H_
#define THIRD_PARTY_NUCLEUS_VENDOR_STATUS_MATCHERS_H_

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace nucleus {

// Matches if a Status or StatusOr object's ok() returns true.
MATCHER(IsOK, "") { return arg.ok(); }

// Matches if a Status or StatusOr object's ok() returns false.
MATCHER(IsNotOK, "") { return !arg.ok(); }

// Matches if a Status or StatusOr object's ok() returns false and that status
// code is expected_code.
MATCHER_P(IsNotOKWithCode, expected_code, "") {
  return (!arg.ok()) &&
         static_cast<int>(arg.code()) == static_cast<int>(expected_code);
}

// Matches if a Status or StatusOr object's ok() returns false and that status
// has an error message containing the string expected_error_message_substring.
MATCHER_P(IsNotOKWithMessage, expected_error_message_substring, "") {
  return (!arg.ok()) &&
         arg.error_message().find(expected_error_message_substring) !=
             std::string::npos;
}

// Matches if a Status or StatusOr object's ok() returns false and that status
// code is expected_code and its error message contains the string
// expected_error_message_substring.
MATCHER_P2(IsNotOKWithCodeAndMessage, expected_code,
           expected_error_message_substring, "") {
  return (!arg.ok()) &&
         static_cast<int>(arg.code()) == static_cast<int>(expected_code) &&
         arg.error_message().find(expected_error_message_substring) !=
             std::string::npos;
}

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_VENDOR_STATUS_MATCHERS_H_
