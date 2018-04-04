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

#ifndef THIRD_PARTY_NUCLEUS_IO_REFERENCE_TEST_H_
#define THIRD_PARTY_NUCLEUS_IO_REFERENCE_TEST_H_

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

using tensorflow::string;

string TestFastaPath() { return GetTestData("test.fasta"); }

typedef std::unique_ptr<GenomeReference> CreateGenomeReferenceFunc(
    const string& fasta_path, int cache_size);

// Tests are parameterized by: reader factory, cache size.
class GenomeReferenceTest : public ::testing::TestWithParam<
    std::pair<CreateGenomeReferenceFunc*, int>> {
 protected:
  void SetUp() override {
    ref_ = (*GetParam().first)(TestFastaPath(), GetParam().second);
  }
  const GenomeReference& Ref() const { return *ref_; }

  std::unique_ptr<const GenomeReference> ref_;
};

typedef GenomeReferenceTest GenomeReferenceDeathTest;

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_REFERENCE_TEST_H_
