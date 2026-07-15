// Phase 5.5e/Path B — Kahan-compensated Conv2D dispatcher.
//
// Wraps a Metal compute pipeline running `conv_kahan_fp32` from
// metal_kernels/conv_kahan_fp32.metal (compiled at runtime via
// `newLibraryWithSource:`). One thread per output element; the
// (kh, kw, c_in) accumulation uses Kahan compensated summation —
// O(ε² · |sum|) per-step error vs O(ε · |sum|) for basic FMA.
// Cross-platform deterministic across reduction orders (Demmel &
// Nguyen ARITH-21 2013, "Fast Reproducible Floating-Point Summation").
//
// Drop-in replacement for `MetalConvSerial::Encode` — same `ConvDesc`
// + buffer layouts. Used to replace MPSGraph conv2D layers where
// non-Kahan reduction drift flips FILTER classes vs Docker (Phase
// 5.5e Path B).
//
// All buffers are FP32 NHWC (input, output) / HWIO (weights).

#pragma once

#include <cstddef>
#include <memory>

#include "deepvariant/native/metal_conv_serial.h"  // reuse ConvDesc

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

class MetalConvKahan {
 public:
  static std::unique_ptr<MetalConvKahan> Create();

  ~MetalConvKahan();

#ifdef __OBJC__
  // Encode one Kahan-compensated Conv2D dispatch into `cmd_buf`.
  // Same parameters and contract as MetalConvSerial::Encode.
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> W, id<MTLBuffer> bias,
              id<MTLBuffer> dst, const ConvDesc& d);
  id<MTLDevice> Device() const;
#endif

  MetalConvKahan(const MetalConvKahan&) = delete;
  MetalConvKahan& operator=(const MetalConvKahan&) = delete;

 private:
  MetalConvKahan();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
