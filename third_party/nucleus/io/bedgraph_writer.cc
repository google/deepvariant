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

#include "third_party/nucleus/io/bedgraph_writer.h"

#include <utility>

#include "absl/log/log.h"
#include "absl/memory/memory.h"
#include "absl/strings/substitute.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

StatusOr<std::unique_ptr<BedGraphWriter>> BedGraphWriter::ToFile(
    const string& bedgraph_path) {
  StatusOr<std::unique_ptr<TextWriter>> text_writer =
      TextWriter::ToFile(bedgraph_path);
  NUCLEUS_RETURN_IF_ERROR(text_writer.status());
  return absl::WrapUnique(new BedGraphWriter(text_writer.ConsumeValueOrDie()));
}

BedGraphWriter::~BedGraphWriter() {
  if (!text_writer_) {
    return;
  }
  ::nucleus::Status status = Close();
  if (!status.ok()) {
    LOG(WARNING) << "Closing BedGraphReader encountered an error";
  }
}

::nucleus::Status BedGraphWriter::Close() {
  if (!text_writer_) {
    return ::nucleus::FailedPrecondition(
        "Cannot close an already closed BedGraphWriter");
  }
  ::nucleus::Status close_status = text_writer_->Close();
  text_writer_ = nullptr;
  return close_status;
}

::nucleus::Status BedGraphWriter::Write(
    const nucleus::genomics::v1::BedGraphRecord& record) {
  if (!text_writer_) {
    return ::nucleus::FailedPrecondition(
        "Cannot write to closed bedgraph stream.");
  }
  NUCLEUS_RETURN_IF_ERROR(text_writer_->Write(
      absl::Substitute("$0\t$1\t$2\t$3\n", record.reference_name(),
                       record.start(), record.end(), record.data_value())));
  return ::nucleus::Status();
}

BedGraphWriter::BedGraphWriter(std::unique_ptr<TextWriter> text_writer)
    : text_writer_(std::move(text_writer)) {}

}  // namespace nucleus
