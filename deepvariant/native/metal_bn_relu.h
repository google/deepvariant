// Phase 5.5f — separate BatchNorm+ReLU dispatcher matching TF/oneDNN's
// non-folded conv→BN→ReLU sequence. Used downstream of MetalConvSerial
// (with folded ReLU disabled) to avoid the FoldConvBn FP32 drift that
// causes ~0.08 % FILTER mismatches on full chr20 vs Docker.
//
// Day-1 PoC measurement: folded conv+BN+ReLU drift up to 93 ULP per
// element on stem_s1a; switching to per-thread c_in-serial FMA conv +
// this kernel reduces max delta to ±2 ULP (76 % bit-exact).

#pragma once

#include <cstddef>
#include <memory>

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandQueue;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct BnReluDesc {
  int B;            // batch size
  int H;            // spatial height
  int W;            // spatial width
  int C;            // channels (mean/var/beta sized to this)
  float eps = 1.0e-3f;  // Keras BN default
  bool relu = true;     // apply ReLU after BN
};

class MetalBnRelu {
 public:
  static std::unique_ptr<MetalBnRelu> Create();
  ~MetalBnRelu();

#ifdef __OBJC__
  // Encode one BN+ReLU dispatch onto `cmd_buf`. Buffers must be FP32
  // on the same device. Sizes:
  //   src   : B * H * W * C  (NHWC, output of preceding raw conv)
  //   mean  : C
  //   var   : C
  //   beta  : C
  //   dst   : B * H * W * C  (NHWC, may alias src for in-place)
  bool Encode(id<MTLCommandBuffer> cmd_buf,
              id<MTLBuffer> src, id<MTLBuffer> mean, id<MTLBuffer> var,
              id<MTLBuffer> beta, id<MTLBuffer> dst,
              const BnReluDesc& d);
  id<MTLDevice> Device() const;
#endif

  MetalBnRelu(const MetalBnRelu&) = delete;
  MetalBnRelu& operator=(const MetalBnRelu&) = delete;

 private:
  MetalBnRelu();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
