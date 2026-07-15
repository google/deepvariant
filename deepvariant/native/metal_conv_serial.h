// Phase 5.5c — deterministic-reduction-order Conv2D dispatcher.
//
// Wraps a Metal compute pipeline running `conv_serial_fp32` from
// metal_kernels/conv_serial_fp32.metal (compiled at runtime via
// `newLibraryWithSource:`). One thread per output element; the
// (kh, kw, c_in) accumulation is sequential FP32 with IEEE FMA —
// bit-identical to TF Eigen's FMA path on x86 AVX-512.
//
// Used to selectively replace MPSGraph `convolution2DWithSourceTensor`
// for layers where MPSGraph's parallel-reduction-order produces
// FILTER-flipping drift vs Docker. See PORT_LOG.md Phase 5.5c.
//
// All buffers are FP32 NHWC (input, output) / HWIO (weights) — matching
// the existing metal_inference.mm conventions.

#pragma once

#include <cstddef>
#include <memory>
#include <string>

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandQueue;
@protocol MTLCommandBuffer;
@protocol MTLComputeCommandEncoder;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct ConvDesc {
  int B;
  int H_in;
  int W_in;
  int C_in;
  int H_out;
  int W_out;
  int C_out;
  int Kh;
  int Kw;
  int stride_h = 1;
  int stride_w = 1;
  int pad_h = 0;        // top zero-pad rows; explicit pad model
  int pad_w = 0;        // left zero-pad cols
  bool relu = true;
};

class MetalConvSerial {
 public:
  // Loads + compiles the kernel against the given device. Returns
  // nullptr on compile or pipeline-state error.
  static std::unique_ptr<MetalConvSerial> Create();

  ~MetalConvSerial();

  // Encode one Conv2D dispatch into `cmd_buf`. Buffers must be valid
  // FP32 (no offset, contiguous) on the same device. Sizes:
  //   src   : B * H_in  * W_in  * C_in   floats
  //   W     : Kh * Kw * C_in * C_out     floats   (HWIO)
  //   bias  : C_out                       floats
  //   dst   : B * H_out * W_out * C_out  floats
  //
  // The encoder must be of compute type and is left in an open state
  // (the caller may queue additional work). Pass nullptr to use a
  // fresh encoder per call (the implementation creates and ends one).
#ifdef __OBJC__
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> W, id<MTLBuffer> bias,
              id<MTLBuffer> dst, const ConvDesc& d);
  id<MTLDevice> Device() const;
#endif

  MetalConvSerial(const MetalConvSerial&) = delete;
  MetalConvSerial& operator=(const MetalConvSerial&) = delete;

 private:
  MetalConvSerial();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// MaxPool 2D dispatcher (Metal compute, NHWC). Used to bridge between
// deterministic Conv2D layers in the stem chain. Max is associative
// in FP32 → output is bit-identical to MPSGraph's maxpool.
struct MaxPoolDesc {
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
};

class MetalMaxPool {
 public:
  static std::unique_ptr<MetalMaxPool> Create();
  ~MetalMaxPool();

#ifdef __OBJC__
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> dst,
              const MaxPoolDesc& d);
#endif

  MetalMaxPool(const MetalMaxPool&) = delete;
  MetalMaxPool& operator=(const MetalMaxPool&) = delete;

 private:
  MetalMaxPool();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
