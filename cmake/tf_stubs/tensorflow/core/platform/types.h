// Minimal TF type stubs — no TF runtime, just aliases for compilation.
#pragma once
#include <cstdint>
#include <string>
namespace tensorflow {
using uint64 = ::uint64_t;
using int64  = ::int64_t;
using uint32 = ::uint32_t;
using string = ::std::string;
}  // namespace tensorflow
