// Phase 5.5e — channel-axis concat dispatcher impl.

#include "deepvariant/native/metal_concat.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr const char* kConcatChannelsFp32Source = R"DVMSL(
#include <metal_stdlib>
using namespace metal;

struct ConcatParams {
    int B;
    int H;
    int W;
    int n_branches;
    int c_size_0;
    int c_size_1;
    int c_size_2;
    int c_size_3;
    int c_total;
};

kernel void concat_channels_fp32(
    constant ConcatParams& P [[ buffer(0) ]],
    device   const float*  src0 [[ buffer(1) ]],
    device   const float*  src1 [[ buffer(2) ]],
    device   const float*  src2 [[ buffer(3) ]],
    device   const float*  src3 [[ buffer(4) ]],
    device   float*        dst  [[ buffer(5) ]],
    uint3 gid                 [[ thread_position_in_grid ]])
{
    const int c_out = (int)gid.x;
    const int hw    = (int)gid.y;
    const int n     = (int)gid.z;
    if (n >= P.B || c_out >= P.c_total || hw >= P.H * P.W) return;
    const int h = hw / P.W;
    const int w = hw % P.W;

    int b = 0;
    int c_local = c_out;
    int c_size = P.c_size_0;
    if (c_local < c_size) {
        b = 0;
    } else {
        c_local -= c_size;
        c_size = P.c_size_1;
        if (c_local < c_size) {
            b = 1;
        } else {
            c_local -= c_size;
            c_size = P.c_size_2;
            if (c_local < c_size) {
                b = 2;
            } else {
                c_local -= c_size;
                b = 3;
            }
        }
    }

    float v;
    const int hw_off = (n * P.H + h) * P.W + w;
    switch (b) {
        case 0: v = src0[hw_off * P.c_size_0 + c_local]; break;
        case 1: v = src1[hw_off * P.c_size_1 + c_local]; break;
        case 2: v = src2[hw_off * P.c_size_2 + c_local]; break;
        default: v = src3[hw_off * P.c_size_3 + c_local]; break;
    }
    dst[hw_off * P.c_total + c_out] = v;
}
)DVMSL";

struct alignas(16) ConcatParamsGpu {
  int B, H, W, n_branches;
  int c_size_0, c_size_1, c_size_2, c_size_3;
  int c_total;
};

}  // namespace

struct MetalConcat::Impl {
  id<MTLDevice> device = nil;
  id<MTLLibrary> library = nil;
  id<MTLFunction> function = nil;
  id<MTLComputePipelineState> pso = nil;
  // Reusable zero buffer for unused branches (concat with n_branches < 4).
  id<MTLBuffer> zero_buffer = nil;
};

MetalConcat::MetalConcat() = default;
MetalConcat::~MetalConcat() = default;

std::unique_ptr<MetalConcat> MetalConcat::Create() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) return nullptr;

    MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
    opts.fastMathEnabled = NO;
    opts.languageVersion = MTLLanguageVersion3_0;

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kConcatChannelsFp32Source];
    id<MTLLibrary> library =
        [device newLibraryWithSource:src options:opts error:&err];
    if (!library) {
      LOG(ERROR) << "MetalConcat::Create: kernel compile failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }
    id<MTLFunction> function =
        [library newFunctionWithName:@"concat_channels_fp32"];
    if (!function) return nullptr;
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:function error:&err];
    if (!pso) {
      LOG(ERROR) << "MetalConcat::Create: PSO failed: "
                 << (err ? err.localizedDescription.UTF8String : "?");
      return nullptr;
    }

    auto self = std::unique_ptr<MetalConcat>(new MetalConcat());
    self->impl_ = std::make_unique<Impl>();
    self->impl_->device = device;
    self->impl_->library = library;
    self->impl_->function = function;
    self->impl_->pso = pso;
    // Allocate a 16 B zero placeholder for unused src buffers.
    self->impl_->zero_buffer = [device newBufferWithLength:16
                                                   options:MTLResourceStorageModeShared];
    if (!self->impl_->zero_buffer) return nullptr;
    memset([self->impl_->zero_buffer contents], 0, 16);
    return self;
  }
}

bool MetalConcat::Encode(id<MTLCommandBuffer> cmd_buf,
                          id<MTLBuffer> src0, id<MTLBuffer> src1,
                          id<MTLBuffer> src2, id<MTLBuffer> src3,
                          id<MTLBuffer> dst, const ConcatDesc& d) {
  if (!cmd_buf || !dst) return false;
  if (d.n_branches < 1 || d.n_branches > 4) return false;

  ConcatParamsGpu params{};
  params.B = d.B;
  params.H = d.H;
  params.W = d.W;
  params.n_branches = d.n_branches;
  params.c_size_0 = d.c_size[0];
  params.c_size_1 = d.n_branches >= 2 ? d.c_size[1] : 0;
  params.c_size_2 = d.n_branches >= 3 ? d.c_size[2] : 0;
  params.c_size_3 = d.n_branches >= 4 ? d.c_size[3] : 0;
  params.c_total = params.c_size_0 + params.c_size_1 +
                    params.c_size_2 + params.c_size_3;

  // Substitute zero buffer for unused inputs (Metal requires non-nil).
  id<MTLBuffer> z = impl_->zero_buffer;
  if (!src0) src0 = z;
  if (!src1) src1 = z;
  if (!src2) src2 = z;
  if (!src3) src3 = z;

  id<MTLComputeCommandEncoder> enc = [cmd_buf computeCommandEncoder];
  [enc setComputePipelineState:impl_->pso];
  [enc setBytes:&params length:sizeof(params) atIndex:0];
  [enc setBuffer:src0 offset:0 atIndex:1];
  [enc setBuffer:src1 offset:0 atIndex:2];
  [enc setBuffer:src2 offset:0 atIndex:3];
  [enc setBuffer:src3 offset:0 atIndex:4];
  [enc setBuffer:dst offset:0 atIndex:5];

  MTLSize grid = MTLSizeMake(static_cast<NSUInteger>(params.c_total),
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

}  // namespace deepvariant
