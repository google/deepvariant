/*
 * Copyright 2024 Google LLC.
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
#include "third_party/nucleus/io/example_writer.h"
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include "absl/base/optimization.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_replace.h"
#include "absl/strings/string_view.h"
#include "third_party/bagz/bag_writer.h"
#include "third_party/re2/re2.h"
#include "tensorflow/core/lib/io/record_writer.h"


namespace nucleus {

ExampleFormat AutodetectFormat(absl::string_view path, ExampleFormat format) {
  if (format != ExampleFormat::kAuto) return format;
  absl::string_view extension;
  RE2::PartialMatch(absl::AsciiStrToLower(path),
                    R"(\.(bagz|tfrecord\.gz)$)",  // extension
                    &extension);
  if (extension == "tfrecord.gz") {
    return ExampleFormat::kTfRecord;
  } else {
    LOG(FATAL) << absl::StrFormat("Unsupported file extension: %s.\n",
                                  extension);
  }
}


class ExampleWriter::Impl {
 public:
  virtual ~Impl() = default;
  virtual bool Add(absl::string_view value, absl::string_view key) = 0;
  virtual bool Close() = 0;
  absl::Status status() { return status_; }

 protected:
  void UpdateStatus(absl::Status status) { status_.Update(std::move(status)); }

 private:
  absl::Status status_;
};



class ExampleWriter::TfRecordImpl : public ExampleWriter::Impl {
 public:
  explicit TfRecordImpl(absl::string_view path) {
    UpdateStatus(tensorflow::Env::Default()->NewWritableFile(std::string(path),
                                                             &tf_file_));
    if (ABSL_PREDICT_FALSE(!status().ok())) {
      LOG(FATAL) << "Failed to create file " << path << "\n ";
      return;
    }

    const tensorflow::io::RecordWriterOptions& options =
        tensorflow::io::RecordWriterOptions::CreateRecordWriterOptions(
            "GZIP");

    tf_ = std::make_unique<tensorflow::io::RecordWriter>(
            tf_file_.get(),
            options);

  }

  bool Add(absl::string_view value, absl::string_view key) override {
    UpdateStatus(tf_->WriteRecord(value));
    return status().ok();
  }

  bool Close() override {
    UpdateStatus(tf_->Close());
    UpdateStatus(tf_file_->Close());
    return status().ok();
  }

 private:
  std::unique_ptr<tensorflow::WritableFile> tf_file_;
  std::unique_ptr<tensorflow::io::RecordWriter> tf_;
};


ExampleWriter::ExampleWriter(absl::string_view path, ExampleFormat format) {
  }
  status_ = impl_->status();
}

ExampleWriter::~ExampleWriter() { Close(); }

bool ExampleWriter::Add(absl::string_view value, absl::string_view key) {
  if (ABSL_PREDICT_FALSE(impl_ == nullptr)) return false;
  if (ABSL_PREDICT_FALSE(!impl_->Add(value, key))) {
    status_.Update(impl_->status());
    return false;
  }
  return true;
}

bool ExampleWriter::Close() {
  if (ABSL_PREDICT_FALSE(impl_ == nullptr)) return false;
  const bool ok = impl_->Close();
  status_.Update(impl_->status());
  impl_.reset();
  return ok;
}

}  // namespace nucleus
