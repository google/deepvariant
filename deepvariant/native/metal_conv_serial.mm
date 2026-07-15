// Phase 5.5c — deterministic-reduction-order Conv2D dispatcher impl.

#include "deepvariant/native/metal_conv_serial.h"

#include <cstdlib>
#include <mutex>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

#include "deepvariant/native/metal_conv_kahan.h"  // Path B: Kahan delegation

namespace deepvariant {

namespace {

// Embedded `metal_kernels/conv_serial_fp32.metal` source. Kept inline so
// the binary is self-contained — the `.metal` file in the source tree
// is the canonical copy, this string is updated by hand to mirror it
// (file is short and stable; PORT_LOG flags any divergence).
constexpr const char* kConvSerialFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct ConvParams {
    int B;
    int H_in;
    int W_in;
    int C_in;
    int H_out;
    int W_out;
    int C_out;
    int Kh;
    int Kw;
    int stride_h;
    int stride_w;
    int pad_h;
    int pad_w;
    int relu;
};

kernel void conv_serial_fp32(
    constant ConvParams& P     [[ buffer(0) ]],
    device   const float* src  [[ buffer(1) ]],
    device   const float* W    [[ buffer(2) ]],
    device   const float* bias [[ buffer(3) ]],
    device   float*       dst  [[ buffer(4) ]],
    uint3 gid                  [[ thread_position_in_grid ]])
{
    const int c_out = (int)gid.x;
    const int hw    = (int)gid.y;
    const int n     = (int)gid.z;
    if (n >= P.B || c_out >= P.C_out || hw >= P.H_out * P.W_out) return;
    const int h_out = hw / P.W_out;
    const int w_out = hw % P.W_out;

    const int h_base = h_out * P.stride_h - P.pad_h;
    const int w_base = w_out * P.stride_w - P.pad_w;

    float acc = 0.0f;
    for (int kh = 0; kh < P.Kh; ++kh) {
        const int h_in = h_base + kh;
        if (h_in < 0 || h_in >= P.H_in) continue;
        for (int kw = 0; kw < P.Kw; ++kw) {
            const int w_in = w_base + kw;
            if (w_in < 0 || w_in >= P.W_in) continue;
            for (int c_in = 0; c_in < P.C_in; ++c_in) {
                const float x = src[
                    ((n * P.H_in + h_in) * P.W_in + w_in) * P.C_in + c_in];
                const float w = W[
                    ((kh * P.Kw + kw) * P.C_in + c_in) * P.C_out + c_out];
                acc = metal::precise::fma(x, w, acc);
            }
        }
    }

    acc += bias[c_out];
    if (P.relu != 0) acc = max(acc, 0.0f);

    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C_out + c_out] = acc;
}
)DVMSL";

// Buffer-0 layout matches `ConvParams` in the kernel above. Keep in
// sync.
struct alignas(16) ConvParamsGpu {
  int B, H_in, W_in, C_in;
  int H_out, W_out, C_out;
  int Kh, Kw;
  int stride_h, stride_w, pad_h, pad_w;
  int relu;
};

}  // namespace

struct MetalConvSerial::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalConvSerial::MetalConvSerial() = default;
MetalConvSerial::~MetalConvSerial() = default;

std::unique_ptr<MetalConvSerial> MetalConvSerial::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
      LOG(ERROR) << "MetalConvSerial::Create: no Metal device";
      return nullptr;
    }
    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) {
      LOG(ERROR) << "MetalConvSerial::Create: cannot create queue";
      return nullptr;
    }

    // Compile the kernel from source. Disable fast-math + FMA
    // contraction so the compiler doesn't try to "optimise" the
    // sequential accumulator into a parallel reduction or rewrite the
    // explicit `metal::precise::fma` calls.
    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kConvSerialFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalConvSerial::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"conv_serial_fp32"];
    if (!function) {
      LOG(ERROR) << "MetalConvSerial::Create: kernel function not found";
      return nullptr;
    }
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalConvSerial::Create: PSO create failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalConvSerial>(new MetalConvSerial());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->queue = queue;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

// Path B (2026-05-10): when DV_METAL_KAHAN=1 is set, delegate ALL
// MetalConvSerial::Encode calls to a singleton MetalConvKahan instance.
// MetalConvKahan implements the SAME ConvDesc + buffer contract as
// MetalConvSerial but accumulates with Kahan-Babuška compensated
// summation (per-thread, sequential), achieving O(ε²·|sum|) reduction
// error vs O(ε·|sum|) for basic FMA. This brings our reduction
// numerically closer to Eigen-x86's chunked-FMA path that Docker uses,
// with the goal of eliminating the residual ~0.02 % FP32 drift at the
// GQ=20 boundary that flips FILTER classes.
//
// Singleton pattern: lazy-init on first call (after env-var check),
// shared across all dispatch sites in the inference path. No API
// changes anywhere — the swap is transparent.
namespace {
std::once_flag g_kahan_init;
std::unique_ptr<MetalConvKahan> g_kahan;
bool g_kahan_enabled = false;

bool KahanEnabled() {
  std::call_once(g_kahan_init, []() {
    const char* env = std::getenv("DV_METAL_KAHAN");
    if (env && env[0] == '1') {
      auto k = MetalConvKahan::Create();
      if (k) {
        g_kahan = std::move(k);
        g_kahan_enabled = true;
        LOG(INFO) << "MetalConvSerial: DV_METAL_KAHAN=1 — delegating all "
                     "Conv2D to MetalConvKahan (compensated summation)";
      } else {
        LOG(WARNING) << "DV_METAL_KAHAN=1 set but MetalConvKahan::Create "
                        "failed — falling back to basic serial FMA";
      }
    }
  });
  return g_kahan_enabled;
}
}  // namespace

bool MetalConvSerial::Encode(id<MTLCommandBuffer> cmd_buf,
                              id<MTLBuffer> src, id<MTLBuffer> W,
                              id<MTLBuffer> bias, id<MTLBuffer> dst,
                              const ConvDesc& d) {
  if (!cmd_buf || !src || !W || !bias || !dst) {
    LOG(ERROR) << "MetalConvSerial::Encode: nil buffer";
    return false;
  }

  // Path B delegation — transparent to all call sites.
  if (KahanEnabled()) {
    return g_kahan->Encode(cmd_buf, src, W, bias, dst, d);
  }

  ConvParamsGpu params{};
  params.B = d.B;
  params.H_in = d.H_in;
  params.W_in = d.W_in;
  params.C_in = d.C_in;
  params.H_out = d.H_out;
  params.W_out = d.W_out;
  params.C_out = d.C_out;
  params.Kh = d.Kh;
  params.Kw = d.Kw;
  params.stride_h = d.stride_h;
  params.stride_w = d.stride_w;
  params.pad_h = d.pad_h;
  params.pad_w = d.pad_w;
  params.relu = d.relu ? 1 : 0;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  [enc setBytes:&params length:sizeof(params) atIndex:0];
  [enc setBuffer:src offset:0 atIndex:1];
  [enc setBuffer:W offset:0 atIndex:2];
  [enc setBuffer:bias offset:0 atIndex:3];
  [enc setBuffer:dst offset:0 atIndex:4];

  // Grid layout: (C_out, H_out * W_out, B). Threadgroup chosen by Metal
  // — `dispatchThreads` automatically clamps to PSO max threads per
  // group and emits boundary checks via the if-guards inside the
  // kernel.
  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(d.C_out),
                              static_cast<NSUInteger>(d.H_out * d.W_out),
                              static_cast<NSUInteger>(d.B));
  NSUInteger w = impl_->pso.threadExecutionWidth;       // SIMD width (≈ 32)
  NSUInteger h = impl_->pso.maxTotalThreadsPerThreadgroup / w;
  if (h == 0) h = 1;
  MTLSize tg = MTLSizeMake(w, h, 1);
  [enc dispatchThreads:grid threadsPerThreadgroup:tg];
  [enc endEncoding];
  return true;
}

id<MTLDevice> MetalConvSerial::Device() const {
  return impl_ ? impl_->device : nil;
}

// ---------------------------------------------------------------------------
// MetalMaxPool — 2-D max pool, NHWC, FP32. Output is bit-identical to
// MPSGraph maxpool (max is associative in FP32; reduction order
// doesn't matter).
// ---------------------------------------------------------------------------

namespace {

constexpr const char* kMaxPoolFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct MaxPoolParams {
    int B;
    int H_in;
    int W_in;
    int C;
    int H_out;
    int W_out;
    int Kh;
    int Kw;
    int stride_h;
    int stride_w;
    int pad_h;
    int pad_w;
};

kernel void maxpool2d_fp32(
    constant MaxPoolParams& P [[ buffer(0) ]],
    device   const float*   src [[ buffer(1) ]],
    device   float*         dst [[ buffer(2) ]],
    uint3 gid                [[ thread_position_in_grid ]])
{
    const int c  = (int)gid.x;
    const int hw = (int)gid.y;
    const int n  = (int)gid.z;
    if (n >= P.B || c >= P.C || hw >= P.H_out * P.W_out) return;
    const int h_out = hw / P.W_out;
    const int w_out = hw % P.W_out;

    const int h_base = h_out * P.stride_h - P.pad_h;
    const int w_base = w_out * P.stride_w - P.pad_w;

    float m = -INFINITY;
    for (int kh = 0; kh < P.Kh; ++kh) {
        const int h_in = h_base + kh;
        if (h_in < 0 || h_in >= P.H_in) continue;
        for (int kw = 0; kw < P.Kw; ++kw) {
            const int w_in = w_base + kw;
            if (w_in < 0 || w_in >= P.W_in) continue;
            const float v = src[
                ((n * P.H_in + h_in) * P.W_in + w_in) * P.C + c];
            if (v > m) m = v;
        }
    }
    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C + c] = m;
}
)DVMSL";

struct alignas(16) MaxPoolParamsGpu {
  int B, H_in, W_in, C, H_out, W_out;
  int Kh, Kw, stride_h, stride_w, pad_h, pad_w;
};

}  // namespace

struct MetalMaxPool::Impl {
  id<MTLDevice> device = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalMaxPool::MetalMaxPool() = default;
MetalMaxPool::~MetalMaxPool() = default;

std::unique_ptr<MetalMaxPool> MetalMaxPool::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) return nullptr;

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kMaxPoolFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalMaxPool::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"maxpool2d_fp32"];
    if (!function) {
      LOG(ERROR) << "MetalMaxPool::Create: function not found";
      return nullptr;
    }
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalMaxPool::Create: PSO create failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalMaxPool>(new MetalMaxPool());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

bool MetalMaxPool::Encode(id<MTLCommandBuffer> cmd_buf,
                           id<MTLBuffer> src, id<MTLBuffer> dst,
                           const MaxPoolDesc& d) {
  if (!cmd_buf || !src || !dst) {
    LOG(ERROR) << "MetalMaxPool::Encode: nil buffer";
    return false;
  }
  MaxPoolParamsGpu p{};
  p.B = d.B;
  p.H_in = d.H_in;
  p.W_in = d.W_in;
  p.C = d.C;
  p.H_out = d.H_out;
  p.W_out = d.W_out;
  p.Kh = d.Kh;
  p.Kw = d.Kw;
  p.stride_h = d.stride_h;
  p.stride_w = d.stride_w;
  p.pad_h = d.pad_h;
  p.pad_w = d.pad_w;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  [enc setBytes:&p length:sizeof(p) atIndex:0];
  [enc setBuffer:src offset:0 atIndex:1];
  [enc setBuffer:dst offset:0 atIndex:2];

  MTLSize grid = MTLSizeMake((NSUInteger)d.C,
                              (NSUInteger)(d.H_out * d.W_out),
                              (NSUInteger)d.B);
  NSUInteger w = impl_->pso.threadExecutionWidth;
  NSUInteger h = impl_->pso.maxTotalThreadsPerThreadgroup / w;
  if (h == 0) h = 1;
  MTLSize tg = MTLSizeMake(w, h, 1);
  [enc dispatchThreads:grid threadsPerThreadgroup:tg];
  [enc endEncoding];
  return true;
}

}  // namespace deepvariant
