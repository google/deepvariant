// Stub — implementations in patches/gfile_macos.cc use POSIX directly.
#pragma once
#include <string>
#include "tensorflow/core/platform/tstring.h"
#include "tensorflow/core/platform/types.h"
namespace tensorflow {
// Empty stub; nucleus::ReadableFile / WritableFile are reimplemented in patches.
struct RandomAccessFile { virtual ~RandomAccessFile() = default; };
struct WritableFile     { virtual ~WritableFile() = default; };
}  // namespace tensorflow
