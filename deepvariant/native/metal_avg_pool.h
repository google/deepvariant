// Phase 5.5e — deterministic-reduction-order AvgPool2D dispatcher.
//
// Wraps `avg_pool_serial_fp32` from
// metal_kernels/avg_pool_serial_fp32.metal. One thread per output
// element; the (kh, kw) accumulation is a strict scalar `for` loop.
//
// All buffers are FP32 NHWC.

#pragma once

#include <cstddef>
#include <memory>

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct AvgPoolDesc {
  int B;
  int H_in;
  int W_in;
  int C;
  int H_out;
  int W_out;
  int Kh;
  int Kw;
  int stride_h = 1;
  int stride_w = 1;
  int pad_h = 0;
  int pad_w = 0;
  // Inception-v3 uses exclude_padding_from_average=True (Keras default
  // for AveragePooling2D with padding='same'). Set to false for the
  // alternative include-padding semantics.
  bool exclude_pad = true;
};

class MetalAvgPool {
 public:
  static std::unique_ptr<MetalAvgPool> Create();
  ~MetalAvgPool();

#ifdef __OBJC__
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> dst,
              const AvgPoolDesc& d);
#endif

  MetalAvgPool(const MetalAvgPool&) = delete;
  MetalAvgPool& operator=(const MetalAvgPool&) = delete;

 private:
  MetalAvgPool();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
