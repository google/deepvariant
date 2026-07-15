// TF-free replacement for third_party/nucleus/io/example_writer.cc.
// Uses our native dv_tfrecord (uncompressed TFRecord) instead of the TF
// io::RecordWriter (GZIP). The native call_variants reader handles
// uncompressed input directly.

#include "third_party/nucleus/io/example_writer.h"

#include <filesystem>
#include <memory>
#include <string>
#include <system_error>

#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "deepvariant/native/tfrecord.h"

namespace nucleus {

// The upstream header forward-declares Impl as a private nested class and
// stores a unique_ptr<Impl>. Provide a concrete definition here so the
// destructor in this translation unit can compile.
class ExampleWriter::Impl {
 public:
  std::unique_ptr<deepvariant::TFRecordWriter> writer;
};

ExampleWriter::ExampleWriter(absl::string_view path, ExampleFormat /*format*/) {
  // Ensure parent directory exists.
  std::filesystem::path p(std::string{path});
  std::error_code ec;
  if (!p.parent_path().empty() &&
      !std::filesystem::is_directory(p.parent_path(), ec)) {
    std::filesystem::create_directories(p.parent_path(), ec);
  }

  impl_ = std::make_unique<Impl>();
  impl_->writer = deepvariant::TFRecordWriter::New(std::string{path});
  if (!impl_->writer) {
    status_ = absl::InternalError(
        absl::StrCat("Failed to open TFRecord writer at ", path));
    impl_.reset();
    return;
  }
  status_ = absl::OkStatus();
}

ExampleWriter::~ExampleWriter() { Close(); }

bool ExampleWriter::Add(absl::string_view value,
                         absl::string_view /*chrom*/, int64_t /*pos*/) {
  if (!impl_ || !impl_->writer) return false;
  // Pass the string_view directly; TFRecordWriter copies into its own
  // coalescing buffer, so materializing an intermediate std::string here
  // would be a redundant copy.
  if (!impl_->writer->WriteRecord(std::string_view(value.data(), value.size()))) {
    status_.Update(absl::InternalError("TFRecord write failed"));
    return false;
  }
  return true;
}

bool ExampleWriter::Close() {
  if (!impl_) return false;
  bool ok = true;
  if (impl_->writer) {
    ok = impl_->writer->Close();
    if (!ok) status_.Update(absl::InternalError("TFRecord close failed"));
  }
  impl_.reset();
  return ok;
}

}  // namespace nucleus
