/*
 * Copyright 2019 Google LLC.
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

#include <memory>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/io/vcf_writer.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/core/status_matchers.h"

namespace nucleus {

namespace {
using genomics::v1::Variant;
}  // namespace

// Read from a vcf file and write to a bcf file. The bcf file should have the
// same contents as the vcf file.
TEST(VcfRoundtripTest, ReadVCFWriteBCF) {
  string input_file = GetTestData("test_sites.vcf");
  string output_file = MakeTempFile("output.bcf.gz");
  auto reader = std::move(
      VcfReader::FromFile(input_file, genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  auto writer = std::move(VcfWriter::ToFile(output_file, reader->Header(),
                                            genomics::v1::VcfWriterOptions())
                              .ValueOrDie());
  std::vector<Variant> expected_variants = as_vector(reader->Iterate());
  for (const auto& v : expected_variants) {
    ASSERT_THAT(writer->Write(v), IsOK());
  }
  writer = nullptr;
  auto output_reader = std::move(
      VcfReader::FromFile(output_file, genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  EXPECT_THAT(as_vector(output_reader->Iterate()),
              testing::Pointwise(EqualsProto(), expected_variants));
}

}  // namespace nucleus
