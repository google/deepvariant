// Phase 5.5e — deterministic global-avg-pool dispatcher impl.

#include "deepvariant/native/metal_global_avg_pool.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr const char* kGlobalAvgPoolFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct GlobalAvgPoolParams {
    int B;
    int H_in;
    int W_in;
    int C;
};

kernel void global_avg_pool_fp32(
    constant GlobalAvgPoolParams& P [[ buffer(0) ]],
    device   const float*         src [[ buffer(1) ]],
    device   float*               dst [[ buffer(2) ]],
    uint2 gid                       [[ thread_position_in_grid ]])
{
    const int c = (int)gid.x;
    const int n = (int)gid.y;
    if (n >= P.B || c >= P.C) return;

    float acc = 0.0f;
    for (int h = 0; h < P.H_in; ++h) {
        for (int w = 0; w < P.W_in; ++w) {
            acc += src[((n * P.H_in + h) * P.W_in + w) * P.C + c];
        }
    }
    const int n_elems = P.H_in * P.W_in;
    dst[n * P.C + c] = acc / (float)n_elems;
}
)DVMSL";

struct alignas(16) GlobalAvgPoolParamsGpu {
  int B, H_in, W_in, C;
};

}  // namespace

struct MetalGlobalAvgPool::Impl {
  id<MTLDevice> device = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalGlobalAvgPool::MetalGlobalAvgPool() = default;
MetalGlobalAvgPool::~MetalGlobalAvgPool() = default;

std::unique_ptr<MetalGlobalAvgPool> MetalGlobalAvgPool::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) return nullptr;

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src =
        [NSString stringWithUTF8String:kGlobalAvgPoolFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalGlobalAvgPool::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"global_avg_pool_fp32"];
    if (!function) return nullptr;
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalGlobalAvgPool::Create: PSO failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalGlobalAvgPool>(
        new MetalGlobalAvgPool());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

bool MetalGlobalAvgPool::Encode(id<MTLCommandBuffer> cmd_buf,
                                 id<MTLBuffer> src, id<MTLBuffer> dst,
                                 const GlobalAvgPoolDesc& d) {
  if (!cmd_buf || !src || !dst) return false;

  GlobalAvgPoolParamsGpu params{};
  params.B = d.B;
  params.H_in = d.H_in;
  params.W_in = d.W_in;
  params.C = d.C;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  [enc setBytes:&params length:sizeof(params) atIndex:0];
  [enc setBuffer:src offset:0 atIndex:1];
  [enc setBuffer:dst offset:0 atIndex:2];

  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(d.C),
                              static_cast<NSUInteger>(d.B),
                              1);
  NSUInteger w = impl_->pso.threadExecutionWidth;
  NSUInteger h = impl_->pso.maxTotalThreadsPerThreadgroup / w;
  if (h == 0) h = 1;
  MTLSize tg = MTLSizeMake(w, h, 1);
  [enc dispatchThreads:grid threadsPerThreadgroup:tg];
  [enc endEncoding];
  return true;
}

}  // namespace deepvariant
