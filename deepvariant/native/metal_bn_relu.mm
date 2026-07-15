// Phase 5.5f — separate BatchNorm+ReLU dispatcher impl.

#include "deepvariant/native/metal_bn_relu.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

// Kept inline so the binary is self-contained. Canonical copy lives at
// metal_kernels/bn_relu_fp32.metal — keep them in sync.
constexpr const char* kBnReluFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct BnReluParams {
    int B;
    int H;
    int W;
    int C;
    float eps;
    int relu;
};

kernel void bn_relu_fp32(
    constant BnReluParams& P     [[ buffer(0) ]],
    device   const float* src    [[ buffer(1) ]],
    device   const float* mean   [[ buffer(2) ]],
    device   const float* var    [[ buffer(3) ]],
    device   const float* beta   [[ buffer(4) ]],
    device   float*       dst    [[ buffer(5) ]],
    uint3 gid                    [[ thread_position_in_grid ]])
{
    const int c = (int)gid.x;
    const int hw = (int)gid.y;
    const int n = (int)gid.z;
    if (n >= P.B || c >= P.C || hw >= P.H * P.W) return;
    const int h = hw / P.W;
    const int w = hw % P.W;
    const int idx = ((n * P.H + h) * P.W + w) * P.C + c;

    const float x = src[idx];
    const float mu = mean[c];
    const float v = var[c];
    const float b = beta[c];

    const float inv_std = 1.0f / metal::precise::sqrt(v + P.eps);
    float y = metal::precise::fma(x - mu, inv_std, b);
    if (P.relu != 0) y = max(y, 0.0f);

    dst[idx] = y;
}
)DVMSL";

struct alignas(16) BnReluParamsGpu {
  int B;
  int H;
  int W;
  int C;
  float eps;
  int relu;
  // Pad to 32 bytes so `setBytes:length:` matches the kernel constant
  // buffer layout.
  int _pad0;
  int _pad1;
};
static_assert(sizeof(BnReluParamsGpu) == 32, "params layout mismatch");

}  // namespace

struct MetalBnRelu::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
};

MetalBnRelu::MetalBnRelu() = default;
MetalBnRelu::~MetalBnRelu() = default;

std::unique_ptr<MetalBnRelu> MetalBnRelu::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
      LOG(ERROR) << "MetalBnRelu::Create: no Metal device";
      return nullptr;
    }
    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) {
      LOG(ERROR) << "MetalBnRelu::Create: cannot create queue";
      return nullptr;
    }

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kBnReluFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalBnRelu::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function = [library newFunctionWithName:@"bn_relu_fp32"];
    if (!function) {
      LOG(ERROR) << "MetalBnRelu::Create: kernel function not found";
      return nullptr;
    }
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalBnRelu::Create: PSO create failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalBnRelu>(new MetalBnRelu());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->queue = queue;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    return self;
  }
}

bool MetalBnRelu::Encode(id<MTLCommandBuffer> cmd_buf,
                          id<MTLBuffer> src, id<MTLBuffer> mean,
                          id<MTLBuffer> var, id<MTLBuffer> beta,
                          id<MTLBuffer> dst, const BnReluDesc& d) {
  if (!cmd_buf || !src || !mean || !var || !beta || !dst) {
    LOG(ERROR) << "MetalBnRelu::Encode: nil buffer";
    return false;
  }

  BnReluParamsGpu params{};
  params.B = d.B;
  params.H = d.H;
  params.W = d.W;
  params.C = d.C;
  params.eps = d.eps;
  params.relu = d.relu ? 1 : 0;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  // Note: kernel reads BnReluParams (24 bytes); pad to 32 bytes for
  // alignment but only the leading bytes are interpreted.
  [enc setBytes:&params length:24 atIndex:0];
  [enc setBuffer:src offset:0 atIndex:1];
  [enc setBuffer:mean offset:0 atIndex:2];
  [enc setBuffer:var offset:0 atIndex:3];
  [enc setBuffer:beta offset:0 atIndex:4];
  [enc setBuffer:dst offset:0 atIndex:5];

  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(d.C),
                              static_cast<NSUInteger>(d.H * d.W),
                              static_cast<NSUInteger>(d.B));
  NSUInteger w = impl_->pso.threadExecutionWidth;
  NSUInteger h = impl_->pso.maxTotalThreadsPerThreadgroup / w;
  if (h == 0) h = 1;
  MTLSize tg = MTLSizeMake(w, h, 1);
  [enc dispatchThreads:grid threadsPerThreadgroup:tg];
  [enc endEncoding];
  return true;
}

id<MTLDevice> MetalBnRelu::Device() const {
  return impl_ ? impl_->device : nil;
}

}  // namespace deepvariant
