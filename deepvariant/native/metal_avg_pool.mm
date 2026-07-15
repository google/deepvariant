// Phase 5.5e — deterministic AvgPool2D dispatcher impl.

#include "deepvariant/native/metal_avg_pool.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

// Embedded `metal_kernels/avg_pool_serial_fp32.metal` source (kept
// inline so the binary is self-contained — the .metal file is the
// canonical copy; updates mirror it).
constexpr const char* kAvgPoolFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct AvgPoolParams {
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
    int exclude_pad;
};

kernel void avgpool2d_fp32(
    constant AvgPoolParams& P [[ buffer(0) ]],
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

    float acc = 0.0f;
    int count = 0;
    for (int kh = 0; kh < P.Kh; ++kh) {
        const int h_in = h_base + kh;
        if (h_in < 0 || h_in >= P.H_in) {
            if (P.exclude_pad == 0) count += P.Kw;
            continue;
        }
        for (int kw = 0; kw < P.Kw; ++kw) {
            const int w_in = w_base + kw;
            if (w_in < 0 || w_in >= P.W_in) {
                if (P.exclude_pad == 0) ++count;
                continue;
            }
            acc += src[
                ((n * P.H_in + h_in) * P.W_in + w_in) * P.C + c];
            ++count;
        }
    }

    const float divisor = (count > 0) ? (float)count : 1.0f;
    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C + c] = acc / divisor;
}
)DVMSL";

struct alignas(16) AvgPoolParamsGpu {
  int B, H_in, W_in, C, H_out, W_out;
  int Kh, Kw, stride_h, stride_w, pad_h, pad_w, exclude_pad;
};

}  // namespace

struct MetalAvgPool::Impl {
  id<MTLDevice> device = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalAvgPool::MetalAvgPool() = default;
MetalAvgPool::~MetalAvgPool() = default;

std::unique_ptr<MetalAvgPool> MetalAvgPool::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) return nullptr;

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kAvgPoolFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalAvgPool::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"avgpool2d_fp32"];
    if (!function) {
      LOG(ERROR) << "MetalAvgPool::Create: function not found";
      return nullptr;
    }
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalAvgPool::Create: PSO create failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalAvgPool>(new MetalAvgPool());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

bool MetalAvgPool::Encode(id<MTLCommandBuffer> cmd_buf,
                           id<MTLBuffer> src, id<MTLBuffer> dst,
                           const AvgPoolDesc& d) {
  if (!cmd_buf || !src || !dst) return false;

  AvgPoolParamsGpu params{};
  params.B = d.B;
  params.H_in = d.H_in;
  params.W_in = d.W_in;
  params.C = d.C;
  params.H_out = d.H_out;
  params.W_out = d.W_out;
  params.Kh = d.Kh;
  params.Kw = d.Kw;
  params.stride_h = d.stride_h;
  params.stride_w = d.stride_w;
  params.pad_h = d.pad_h;
  params.pad_w = d.pad_w;
  params.exclude_pad = d.exclude_pad ? 1 : 0;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  [enc setBytes:&params length:sizeof(params) atIndex:0];
  [enc setBuffer:src offset:0 atIndex:1];
  [enc setBuffer:dst offset:0 atIndex:2];

  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(d.C),
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

}  // namespace deepvariant
