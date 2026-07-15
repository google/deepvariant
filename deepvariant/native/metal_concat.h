// Phase 5.5e — channel-axis concat dispatcher (NHWC FP32).
//
// One thread per output element; pure data movement. Up to 4 input
// branches (matches Inception-v3 max-branch count).

#pragma once

#include <cstddef>
#include <memory>

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct ConcatDesc {
  int B;
  int H;
  int W;
  int n_branches;          // 1..4
  int c_size[4];           // channel count per branch (unused entries 0)
  // c_total computed by Encode().
};

class MetalConcat {
 public:
  static std::unique_ptr<MetalConcat> Create();
  ~MetalConcat();

#ifdef __OBJC__
  // Pass nullptr for unused branches when n_branches < 4.
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src0, id<MTLBuffer> src1,
              id<MTLBuffer> src2, id<MTLBuffer> src3,
              id<MTLBuffer> dst, const ConcatDesc& d);
#endif

  MetalConcat(const MetalConcat&) = delete;
  MetalConcat& operator=(const MetalConcat&) = delete;

 private:
  MetalConcat();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
