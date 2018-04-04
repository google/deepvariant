/*
 * Copyright 2018 Google Inc.
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

#include "htslib/hts.h"

#include "tensorflow/core/platform/test.h"
#include "absl/strings/str_cat.h"
#include "htslib/hts_log.h"
#include "htslib/kstring.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/port.h"

using absl::StrCat;

namespace nucleus {

// This tests minimal htslib functionality.
TEST(FileTest, HtsOpen) {
  hts_set_log_level(HTS_LOG_TRACE);

  string path = GetTestData("test.fasta");
  errno = 0;
  htsFile* f = hts_open(path.c_str(), "rb");
  EXPECT_EQ(errno, 0);
  ASSERT_NE(f, nullptr);

  kstring_t str = {0, 0, nullptr};
  int n = hts_getline(f, '\n', &str);
  EXPECT_GT(n, 0);
  EXPECT_EQ(ks_len(&str), 5);
  EXPECT_EQ(ks_str(&str), string(">chrM"));
  free(str.s);

  int close_ok = hts_close(f);
  EXPECT_EQ(close_ok, 0);
}

}  // namespace nucleus
