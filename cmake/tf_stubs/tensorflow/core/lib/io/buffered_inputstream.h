// Stub — ReadableFile is reimplemented in patches/gfile_macos.cc.
#pragma once
#include <string>
#include "tensorflow/core/platform/file_system.h"
namespace tensorflow { namespace io {
struct RandomAccessInputStream { explicit RandomAccessInputStream(RandomAccessFile*, bool) {} };
struct BufferedInputStream {
  BufferedInputStream(RandomAccessInputStream*, size_t, bool) {}
  bool ReadLine(std::string*) { return false; }
};
}}  // namespace tensorflow::io
