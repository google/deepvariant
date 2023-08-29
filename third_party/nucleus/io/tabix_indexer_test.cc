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

#include "third_party/nucleus/io/tabix_indexer.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/io/vcf_writer.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

using ::testing::Test;

constexpr char kVcfIndexSamplesFilename[] = "test_samples.vcf.gz";

TEST(TabixIndexerTest, IndexBuildsCorrectly) {
  string output_filename = MakeTempFile("test_samples.vcf.gz");
  string output_tabix_index = output_filename + ".tbi";

  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(GetTestData(kVcfIndexSamplesFilename),
                                   nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());

  nucleus::genomics::v1::VcfWriterOptions writer_options;
  std::unique_ptr<VcfWriter> writer =
      std::move(nucleus::VcfWriter::ToFile(output_filename, reader->Header(),
                                           writer_options)
                    .ValueOrDie());

  auto variants = nucleus::as_vector(reader->Iterate());
  for (const auto& v : variants) {
    NUCLEUS_CHECK_OK(writer->Write(v));
  }

  EXPECT_THAT(TbxIndexBuild(output_filename), IsOK());
  EXPECT_THAT(reader->Query(MakeRange("chr3", 14318, 14319)), IsOK());
}

TEST(CSIIndexerTest, IndexBuildsCorrectly) {
  string output_filename = MakeTempFile("test_samples.vcf.gz");
  string output_csi_index = output_filename + ".csi";

  std::unique_ptr<nucleus::VcfReader> reader = std::move(
      nucleus::VcfReader::FromFile(GetTestData(kVcfIndexSamplesFilename),
                                   nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());

  nucleus::genomics::v1::VcfWriterOptions writer_options;
  std::unique_ptr<VcfWriter> writer =
      std::move(nucleus::VcfWriter::ToFile(output_filename, reader->Header(),
                                           writer_options)
                    .ValueOrDie());

  auto variants = nucleus::as_vector(reader->Iterate());
  for (const auto& v : variants) {
    NUCLEUS_CHECK_OK(writer->Write(v));
  }

  EXPECT_THAT(CSIIndexBuild(output_filename, 14), IsOK());
  EXPECT_THAT(tensorflow::Env::Default()->FileExists(output_csi_index), IsOK());
  EXPECT_THAT(reader->Query(MakeRange("chr3", 14318, 14319)), IsOK());
}
}  // namespace nucleus
