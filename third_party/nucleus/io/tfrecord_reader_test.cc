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
#include <string>

#include "third_party/nucleus/io/tfrecord_reader.h"

#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace nucleus {

TEST(TFRecordReaderTest, Simple) {
  std::unique_ptr<TFRecordReader> reader = TFRecordReader::New(
      GetTestData("test_likelihoods.vcf.golden.tfrecord"), "");
  ASSERT_NE(reader, nullptr);

  ASSERT_TRUE(reader->GetNext());

  tensorflow::tstring s = reader->record();

  nucleus::genomics::v1::Variant v;
  v.ParseFromArray(s.data(), s.size());

  ASSERT_EQ("Chr1", v.reference_name());

  reader->Close();
}


TEST(TFRecordReaderTest, NotFound) {
  std::unique_ptr<TFRecordReader> reader =
      TFRecordReader::New(GetTestData("not_found.tfrecord"), "");
  ASSERT_EQ(reader, nullptr);
}

}  // namespace nucleus

