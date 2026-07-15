// Phase 5.5e — deterministic global-avg-pool dispatcher.
//
// Reduces NHWC (B, H_in, W_in, C) → (B, C) by averaging over the
// spatial volume. One thread per output element (n, c). Per-thread
// strict-serial accumulation ensures bit-determinism.

#pragma once

#include <cstddef>
#include <memory>

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct GlobalAvgPoolDesc {
  int B;
  int H_in;
  int W_in;
  int C;
};

class MetalGlobalAvgPool {
 public:
  static std::unique_ptr<MetalGlobalAvgPool> Create();
  ~MetalGlobalAvgPool();

#ifdef __OBJC__
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> dst,
              const GlobalAvgPoolDesc& d);
#endif

  MetalGlobalAvgPool(const MetalGlobalAvgPool&) = delete;
  MetalGlobalAvgPool& operator=(const MetalGlobalAvgPool&) = delete;

 private:
  MetalGlobalAvgPool();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
