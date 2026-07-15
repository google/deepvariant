// Stub — TFRecord writer is reimplemented in patches/tfrecord_writer_macos.cc.
#pragma once
#include "tensorflow/core/platform/file_system.h"
namespace tensorflow { namespace io {
struct RecordWriterOptions {};
struct RecordWriter {
  RecordWriter(WritableFile*, const RecordWriterOptions& = {}) {}
  void Flush() {}
  void Close() {}
};
}}  // namespace tensorflow::io
