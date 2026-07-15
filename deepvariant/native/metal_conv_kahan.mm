// Phase 5.5e/Path B — Kahan-compensated Conv2D dispatcher impl.
//
// Mirrors `metal_conv_serial.mm` (same params, same dispatch shape)
// but with Kahan compensation in the inner accumulation loop. The
// embedded kernel source must match
// `metal_kernels/conv_kahan_fp32.metal` byte-for-byte (the .metal file
// is the canonical copy).

#include "deepvariant/native/metal_conv_kahan.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr const char* kConvKahanFp32Source = R"DVMSL(
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

kernel void conv_kahan_fp32(
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

    float sum = 0.0f;
    float c = 0.0f;
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
                const float y = metal::precise::fma(x, w, -c);
                const float t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
        }
    }

    sum += bias[c_out];
    if (P.relu != 0) sum = max(sum, 0.0f);

    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C_out + c_out] = sum;
}
)DVMSL";

struct alignas(16) ConvParamsGpu {
  int B, H_in, W_in, C_in;
  int H_out, W_out, C_out;
  int Kh, Kw;
  int stride_h, stride_w, pad_h, pad_w;
  int relu;
};

}  // namespace

struct MetalConvKahan::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalConvKahan::MetalConvKahan() = default;
MetalConvKahan::~MetalConvKahan() = default;

std::unique_ptr<MetalConvKahan> MetalConvKahan::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
      LOG(ERROR) << "MetalConvKahan::Create: no Metal device";
      return nullptr;
    }
    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) {
      LOG(ERROR) << "MetalConvKahan::Create: cannot create queue";
      return nullptr;
    }

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kConvKahanFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalConvKahan::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"conv_kahan_fp32"];
    if (!function) {
      LOG(ERROR) << "MetalConvKahan::Create: kernel function not found";
      return nullptr;
    }
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalConvKahan::Create: PSO create failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalConvKahan>(new MetalConvKahan());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->queue = queue;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

bool MetalConvKahan::Encode(id<MTLCommandBuffer> cmd_buf,
                             id<MTLBuffer> src, id<MTLBuffer> W,
                             id<MTLBuffer> bias, id<MTLBuffer> dst,
                             const ConvDesc& d) {
  if (!cmd_buf || !src || !W || !bias || !dst) {
    LOG(ERROR) << "MetalConvKahan::Encode: nil buffer";
    return false;
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

  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(d.C_out),
                              static_cast<NSUInteger>(d.H_out * d.W_out),
                              static_cast<NSUInteger>(d.B));
  NSUInteger w = impl_->pso.threadExecutionWidth;
  NSUInteger h = impl_->pso.maxTotalThreadsPerThreadgroup / w;
  if (h == 0) h = 1;
  MTLSize tg = MTLSizeMake(w, h, 1);
  [enc dispatchThreads:grid threadsPerThreadgroup:tg];
  [enc endEncoding];
  return true;
}

id<MTLDevice> MetalConvKahan::Device() const {
  return impl_ ? impl_->device : nil;
}

}  // namespace deepvariant
