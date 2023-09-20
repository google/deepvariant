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

#include "deepvariant/postprocess_variants.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/io/compression.h"
#include "tensorflow/core/lib/io/record_reader.h"
#include "tensorflow/core/lib/io/record_writer.h"

namespace learning {
namespace genomics {
namespace deepvariant {

namespace {

void SortSingleSiteCalls(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    std::vector<CallVariantsOutput>* calls) {
  std::vector<CallVariantsOutput> output;
  if (calls->empty()) {
    return;
  }
  //   Create the mapping from from contig to pos_in_fasta.
  std::map<std::string, int> contig_name_to_pos_in_fasta =
      nucleus::MapContigNameToPosInFasta(contigs);
  std::stable_sort(calls->begin(), calls->end(),
            [&contig_name_to_pos_in_fasta](const CallVariantsOutput& a,
                                           const CallVariantsOutput& b) {
              return nucleus::CompareVariants(a.variant(), b.variant(),
                                              contig_name_to_pos_in_fasta);
            });
}

}  // namespace

void ProcessSingleSiteCallTfRecords(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<std::string>& tfrecord_paths,
    const string& output_tfrecord_path) {
  std::vector<CallVariantsOutput> single_site_calls;
  tensorflow::Env* env = tensorflow::Env::Default();
  for (const string& tfrecord_path : tfrecord_paths) {
    std::unique_ptr<tensorflow::RandomAccessFile> read_file;
    TF_CHECK_OK(env->NewRandomAccessFile(tfrecord_path, &read_file));
    const char* const option = nucleus::EndsWith(tfrecord_path, ".gz")
                                   ? tensorflow::io::compression::kGzip
                                   : tensorflow::io::compression::kNone;
    tensorflow::io::RecordReader reader(
        read_file.get(),
        tensorflow::io::RecordReaderOptions::CreateRecordReaderOptions(option));

    std::uint64_t offset = 0;
    tensorflow::tstring data;
    LOG(INFO) << "Read from: " << tfrecord_path;
    while (reader.ReadRecord(&offset, &data).ok()) {
      CallVariantsOutput single_site_call;
      QCHECK(single_site_call.ParseFromArray(data.data(), data.length()))
          << "Failed to parse CallVariantsOutput";
      // Here we assume each variant has only 1 call.
      QCHECK_EQ(single_site_call.variant().calls_size(), 1);
      single_site_calls.push_back(single_site_call);
    }
    if (tfrecord_paths.size() > 1) {
      LOG(INFO) << "Done reading: " << tfrecord_path
                << ". #entries in single_site_calls = "
                << std::to_string(single_site_calls.size());
    }
  }
  LOG(INFO) << "Total #entries in single_site_calls = "
            << std::to_string(single_site_calls.size());
  VLOG(3) << "Start SortSingleSiteCalls";
  SortSingleSiteCalls(contigs, &single_site_calls);
  VLOG(3) << "Done SortSingleSiteCalls";

  // Write sorted calls to output_tfrecord_path.
  std::unique_ptr<tensorflow::WritableFile> output_file;
  TF_CHECK_OK(tensorflow::Env::Default()->NewWritableFile(output_tfrecord_path,
                                                          &output_file));
  tensorflow::io::RecordWriter output_writer(output_file.get());
  for (const auto& single_site_call : single_site_calls) {
    tensorflow::Status writer_status =
        output_writer.WriteRecord(single_site_call.SerializeAsString());
    QCHECK(writer_status.ok())
        << "Failed to write serialized proto to output_writer. "
        << "Status = " << writer_status.error_message();
  }
  TF_CHECK_OK(output_writer.Flush()) << "Failed to flush the output writer.";
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
