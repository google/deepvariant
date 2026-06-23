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
#include <filesystem>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include "absl/base/optimization.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/strings/ascii.h"
#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_replace.h"
#include "absl/strings/string_view.h"
#include "re2/re2.h"
#include "tensorflow/core/lib/io/record_writer.h"


namespace nucleus {

// Returns the TfRecord compression codec implied by the file-name suffix:
// "SNAPPY" for a ".snappy" extension, otherwise "GZIP". This mirrors the
// reader-side autodetection in dv_utils.compression_type_for_examples_path so a
// single source of truth (the file name) drives both ends.
std::string CompressionTypeForPath(absl::string_view path) {
  return absl::EndsWith(absl::AsciiStrToLower(path), ".snappy") ? "SNAPPY"
                                                                : "GZIP";
}

ExampleFormat AutodetectFormat(absl::string_view path,
                               ExampleFormat format) {
  if (format != ExampleFormat::kAuto) return format;
  std::string path_str = absl::AsciiStrToLower(path);
  std::string extension;
  // Replace shard strings and the compression extension (.gz or .snappy).
  RE2::GlobalReplace(&path_str,
                     R"(\@[0-9]+|-\*?\d*-of-\*?\d*|\.gz|\.snappy)", "");
  RE2::PartialMatch(path_str,
                    R"(\.(bagz|tfrecords?)$)",  // extension
                    &extension);
  if (extension == "tfrecord" || extension == "tfrecords") {
    return ExampleFormat::kTfRecord;
  } else {
   LOG(FATAL) << absl::StrFormat("Unsupported file extension: %s.\n",
                                 std::filesystem::path(path).filename());
  }
}


class ExampleWriter::Impl {
 public:
  virtual ~Impl() = default;
  virtual bool Add(absl::string_view value,
                   absl::string_view chrom,
                   int64_t pos) = 0;
  virtual bool Close() = 0;
  absl::Status status() { return status_; }

 protected:
  void UpdateStatus(absl::Status status) { status_.Update(std::move(status)); }

 private:
  absl::Status status_;
};



class ExampleWriter::TfRecordImpl : public ExampleWriter::Impl {
 public:
  TfRecordImpl(absl::string_view path, int compression_level) {
    UpdateStatus(tensorflow::Env::Default()->NewWritableFile(std::string(path),
                                                             &tf_file_));
    if (ABSL_PREDICT_FALSE(!status().ok())) {
      LOG(FATAL) << "Failed to create file " << path << "\n ";
      return;
    }

    // The codec is inferred from the file-name suffix (".snappy" -> SNAPPY,
    // otherwise GZIP) so it always agrees with what the readers autodetect.
    const std::string codec = CompressionTypeForPath(path);
    tensorflow::io::RecordWriterOptions options =
        tensorflow::io::RecordWriterOptions::CreateRecordWriterOptions(codec);
    // The level only applies to the zlib-based codecs (GZIP/ZLIB); a negative
    // value leaves the library default (Z_DEFAULT_COMPRESSION) in place. Snappy
    // has no notion of a level, so a level set alongside a ".snappy" path is
    // ignored -- warn rather than silently dropping it.
    if (compression_level >= 0) {
      if (codec == "SNAPPY") {
        LOG(WARNING) << "Ignoring compression level " << compression_level
                     << " for Snappy output " << path
                     << " (Snappy does not support compression levels).";
      } else {
        options.zlib_options.compression_level = compression_level;
      }
    }

    tf_ = std::make_unique<tensorflow::io::RecordWriter>(
            tf_file_.get(),
            options);
  }

  bool Add(absl::string_view value,
           absl::string_view chrom,
           int64_t pos) override {
    // chrom and pos are unused.
    UpdateStatus(tf_->WriteRecord(value));
    return status().ok();
  }

  bool Close() override {
    UpdateStatus(tf_->Close());
    UpdateStatus(tf_file_->Close());

    return status().ok();
  }

 private:
  std::unique_ptr<tensorflow::io::RecordWriter> tf_;
  std::unique_ptr<tensorflow::WritableFile> tf_file_;
};


ExampleWriter::ExampleWriter(absl::string_view path,
                             ExampleFormat format,
                             int compression_level) {
  std::filesystem::path p = std::filesystem::path(path);
  if (!std::filesystem::is_directory(p.parent_path())
      && p.parent_path() != "") {
  status_ = std::filesystem::create_directories(p.parent_path()) ?
     absl::OkStatus()
     : absl::InternalError("Failed to create directories.");
  }
  switch (format = AutodetectFormat(path, format)) {
    case ExampleFormat::kAuto:
      LOG(INFO) << "Autodetect failure.";
      status_ = absl::InternalError("Autodetect failure.");
      return;
    case ExampleFormat::kTfRecord:
      LOG(INFO) << "Writing output using TfRecord";
      impl_ = std::make_unique<TfRecordImpl>(path, compression_level);
      break;
    case ExampleFormat::kBagz:
      break;
  }
  status_ = impl_->status();
}

ExampleWriter::~ExampleWriter() { Close(); }

bool ExampleWriter::Add(absl::string_view value,
                        absl::string_view chrom,
                        int64_t pos) {
  if (ABSL_PREDICT_FALSE(impl_ == nullptr)) return false;
  if (ABSL_PREDICT_FALSE(!impl_->Add(value, chrom, pos))) {
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
