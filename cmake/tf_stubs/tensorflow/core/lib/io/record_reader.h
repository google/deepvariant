// Stub — TFRecord reader is reimplemented in patches/tfrecord_reader_macos.cc.
#pragma once
#include "tensorflow/core/platform/types.h"
#include "tensorflow/core/platform/tstring.h"
namespace tensorflow { namespace io {
struct RecordReaderOptions {
  static RecordReaderOptions CreateRecordReaderOptions(const std::string&) {
    return RecordReaderOptions{};
  }
};
struct RecordReader { RecordReader(void*, const RecordReaderOptions&) {} };
}}  // namespace tensorflow::io
