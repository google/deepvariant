// MPSGraph implementation of DeepVariant Inception-v3, mirroring
// tools/conversion/inception_v3_mil.py.
//
// All conv+BN pairs are fused on CPU at graph-build time:
//   scale[o]  = 1 / sqrt(var[o] + epsilon)         (gamma is frozen at 1)
//   offset[o] = beta[o] - mean[o] * scale[o]
//   W'[o,i,h,w] = W[o,i,h,w] * scale[o]
// then a single Conv2D + bias-add is emitted to MPSGraph.
//
// MPSGraph data layout: NCHW. We transpose (N,100,221,7) → (N,7,100,221)
// at the input. Concat axis is 1 (channels in NCHW).

#include "deepvariant/native/metal_inference.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShadersGraph/MetalPerformanceShadersGraph.h>
#import <MetalPerformanceShadersGraph/MPSGraphImToColOps.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "absl/log/log.h"
#include "absl/strings/str_split.h"
#include "deepvariant/native/dv_weights.h"
#include "deepvariant/native/metal_bn_relu.h"
#include "deepvariant/native/metal_conv_serial.h"
#include "deepvariant/native/metal_avg_pool.h"
#include "deepvariant/native/metal_concat.h"
#include "deepvariant/native/metal_det_mixed.h"
#include "deepvariant/native/metal_global_avg_pool.h"

namespace deepvariant {

namespace {

// Keras `BatchNormalization` defaults to epsilon=1e-3 (NOT 1e-4) and that
// is the value Inception-v3 SavedModels were trained with. Using 1e-4
// produces a subtle scale mismatch on channels where var is small enough
// that the +eps term changes magnitude — large enough to flip the sign
// of post-ReLU activations on those channels, which manifests as
// channel-level mismatch vs TF reference.
constexpr float kBNEpsilon = 1e-3f;

// Build the layer-N variable name used by extract_weights.py.
std::string AttrCpp(int n, const char* attr) {
  return std::string("layer_with_weights-") + std::to_string(n) +
         "/" + attr + "/.ATTRIBUTES/VARIABLE_VALUE";
}

// Fold a Conv (HWIO, FP32) and a BN (gamma=1, beta, mean, var, epsilon)
// into a fused (W', b') pair in HWIO layout (TF-native — no host-side
// transpose; passed straight into MPSGraph with weightsLayout=HWIO).
struct FusedConv {
  std::vector<float> weights_hwio;  // [H, W, I, O]
  std::vector<float> bias;          // [O]
  int O = 0, I = 0, H = 0, W = 0;
};

FusedConv FoldConvBn(const DvwWeights& dvw, int conv_n, int bn_n) {
  const auto* k = dvw.Get(AttrCpp(conv_n, "kernel"));
  const auto* beta = dvw.Get(AttrCpp(bn_n, "beta"));
  const auto* mean = dvw.Get(AttrCpp(bn_n, "moving_mean"));
  const auto* var = dvw.Get(AttrCpp(bn_n, "moving_variance"));
  if (!k || !beta || !mean || !var) {
    LOG(ERROR) << "FoldConvBn(conv=" << conv_n << ", bn=" << bn_n
               << "): missing weight tensor";
    return {};
  }
  if (k->shape.size() != 4u) {
    LOG(ERROR) << "kernel for layer " << conv_n
               << " has rank " << k->shape.size() << " (need 4)";
    return {};
  }
  // Source layout is HWIO (TF Keras convention): shape = (H, W, I, O).
  const int Hk = k->shape[0];
  const int Wk = k->shape[1];
  const int Ik = k->shape[2];
  const int Ok = k->shape[3];
  if (beta->shape.size() != 1u || (int)beta->shape[0] != Ok ||
      mean->shape.size() != 1u || (int)mean->shape[0] != Ok ||
      var->shape.size() != 1u || (int)var->shape[0] != Ok) {
    LOG(ERROR) << "BN params shape mismatch for conv=" << conv_n
               << " bn=" << bn_n;
    return {};
  }
  FusedConv out;
  out.O = Ok;
  out.I = Ik;
  out.H = Hk;
  out.W = Wk;

  // scale[o] = 1 / sqrt(var[o] + eps);  offset[o] = beta[o] - mean[o]*scale[o]
  std::vector<float> scale(Ok), offset(Ok);
  for (int o = 0; o < Ok; ++o) {
    scale[o] = 1.0f / std::sqrt(var->data[o] + kBNEpsilon);
    offset[o] = beta->data[o] - mean->data[o] * scale[o];
  }
  out.bias = std::move(offset);

  // Native HWIO; multiply by scale[o] along the O axis. Optionally
  // flip H and W axes ("true convolution" vs cross-correlation) — see
  // the diagnostic experiment in PORT_LOG. TF uses cross-correlation.
  out.weights_hwio.resize((size_t)Hk * Wk * Ik * Ok);
  const bool flip_spatial = false;  // TF/MPS conv is cross-correlation, no flip
  for (size_t h = 0; h < (size_t)Hk; ++h) {
    const size_t h_src = flip_spatial ? (Hk - 1 - h) : h;
    for (size_t w = 0; w < (size_t)Wk; ++w) {
      const size_t w_src = flip_spatial ? (Wk - 1 - w) : w;
      for (size_t i = 0; i < (size_t)Ik; ++i) {
        for (size_t o = 0; o < (size_t)Ok; ++o) {
          const size_t dst_idx = ((h * Wk + w) * Ik + i) * Ok + o;
          const size_t src_idx = ((h_src * Wk + w_src) * Ik + i) * Ok + o;
          out.weights_hwio[dst_idx] = k->data[src_idx] * scale[o];
        }
      }
    }
  }
  return out;
}

// MPSGraph constant-tensor helper.
//
// IMPORTANT: must use `[[NSData alloc] initWithBytes:length:]` rather
// than `[NSData dataWithBytes:length:]` here. The latter returns an
// AUTORELEASED NSData; when the autoreleasepool from `Create()` drains
// (i.e. before the first `PredictAtTap()` runs), MPSGraph's internal
// reference to the bytes becomes a dangling pointer and the constant
// tensor reads garbage. The +1-retained alloc/init form keeps the
// NSData alive for as long as ARC tracks it through the MPSGraphTensor
// reference graph, surviving past the build-time pool drain.
//
// This is the Phase 5.5a root cause: months of mysterious channel-
// permutation behaviour traced to autoreleased NSData in the conv
// weight constants.
MPSGraphTensor* ConstFloat32(MPSGraph* g, const float* data,
                             NSArray<NSNumber*>* shape, NSString* name) {
  size_t n = 1;
  for (NSNumber* d in shape) n *= [d unsignedIntegerValue];
  NSData* nsdata = [[NSData alloc] initWithBytes:data
                                          length:n * sizeof(float)];
  return [g constantWithData:nsdata
                       shape:shape
                    dataType:MPSDataTypeFloat32];
}

// Build a Conv2D + bias-add via MPSGraph's native
// `convolution2DWithSourceTensor:` (NHWC + HWIO).
//
// Verified bit-exact for the exact stem_s1a shape (input 100×221×7,
// kernel 3×3 stride-2 valid 7→32) by `microtest_metal` (Phase 5.5a
// investigation, Test 6 — known-pattern weights and sparse input,
// hand-computed expected output, max-abs = 0). Earlier reports of a
// channel-permutation bug here were artifacts of an unrelated stale
// shape-mismatch path in `debug_metal`'s TapList — not a real
// MPSGraph issue.
MPSGraphTensor* AddConv(MPSGraph* g, MPSGraphTensor* x,
                        const FusedConv& fc,
                        int stride_y, int stride_x,
                        bool same_padding,  // true = "same", false = "valid"
                        NSString* name) {
  NSArray* w_shape = @[@(fc.H), @(fc.W), @(fc.I), @(fc.O)];
  MPSGraphTensor* W = ConstFloat32(g, fc.weights_hwio.data(), w_shape,
                                   [name stringByAppendingString:@"_w"]);
  MPSGraphTensor* b = ConstFloat32(g, fc.bias.data(), @[@(fc.O)],
                                   [name stringByAppendingString:@"_b"]);

  MPSGraphConvolution2DOpDescriptor* desc =
      [MPSGraphConvolution2DOpDescriptor
          descriptorWithStrideInX:stride_x
                        strideInY:stride_y
                  dilationRateInX:1
                  dilationRateInY:1
                           groups:1
                     paddingStyle:same_padding
                                      ? MPSGraphPaddingStyleTF_SAME
                                      : MPSGraphPaddingStyleTF_VALID
                       dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                    weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
  MPSGraphTensor* y = [g convolution2DWithSourceTensor:x
                                          weightsTensor:W
                                             descriptor:desc
                                                   name:name];
  // Bias broadcast along channel dim. Bias shape (O,) needs reshape to
  // (1, 1, 1, O) for NHWC broadcasting.
  MPSGraphTensor* b_reshaped = [g reshapeTensor:b
                                       withShape:@[@1, @1, @1, @(fc.O)]
                                            name:[name stringByAppendingString:@"_br"]];
  return [g additionWithPrimaryTensor:y
                       secondaryTensor:b_reshaped
                                  name:[name stringByAppendingString:@"_bias"]];
}

// Phase 5.5f: when true, CBR/AvgCBR build conv + primitive-op BN + ReLU
// (raw kernel weights, BN as separate MPSGraph ops). Gated on
// DV_METAL_UNFOLDED_BN environment variable; set once at MetalInception::
// Create().
static bool g_unfold_bn_for_graph = false;

// Phase 5.5f Conv→BN→ReLU using primitive MPSGraph ops (no FoldConvBn).
//   conv_raw      = conv2D(x, raw_kernel)
//   bn(z) = (z - mean) * inv_std + beta,  inv_std = 1 / sqrt(var + eps)
//   y             = relu(bn(conv_raw))
// inv_std is precomputed host-side (same value as TF computes; only the
// reduction-order through the conv differs from the folded path).
MPSGraphTensor* CBRUnfolded(MPSGraph* g, MPSGraphTensor* x,
                              const DvwWeights& dvw,
                              int conv_n, int bn_n,
                              int stride_y, int stride_x,
                              bool same_padding,
                              NSString* name) {
  const auto* k = dvw.Get(AttrCpp(conv_n, "kernel"));
  const auto* beta = dvw.Get(AttrCpp(bn_n, "beta"));
  const auto* mean = dvw.Get(AttrCpp(bn_n, "moving_mean"));
  const auto* var = dvw.Get(AttrCpp(bn_n, "moving_variance"));
  if (!k || !beta || !mean || !var || k->shape.size() != 4u) {
    LOG(ERROR) << "CBRUnfolded: missing weight for conv=" << conv_n
               << " bn=" << bn_n;
    return nullptr;
  }
  const int Hk = k->shape[0], Wk = k->shape[1];
  const int Ik = k->shape[2], Ok = k->shape[3];

  // Raw conv2D, no bias.
  NSArray* w_shape = @[@(Hk), @(Wk), @(Ik), @(Ok)];
  MPSGraphTensor* W = ConstFloat32(g, k->data, w_shape,
                                   [name stringByAppendingString:@"_w"]);
  MPSGraphConvolution2DOpDescriptor* desc =
      [MPSGraphConvolution2DOpDescriptor
          descriptorWithStrideInX:stride_x
                        strideInY:stride_y
                  dilationRateInX:1
                  dilationRateInY:1
                           groups:1
                     paddingStyle:same_padding ? MPSGraphPaddingStyleTF_SAME
                                               : MPSGraphPaddingStyleTF_VALID
                       dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                    weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
  MPSGraphTensor* conv = [g convolution2DWithSourceTensor:x
                                             weightsTensor:W
                                                descriptor:desc
                                                      name:name];

  // Primitive BN: (conv - mean) * inv_std + beta, with inv_std precomputed
  // host-side. Each tensor is a (1, 1, 1, Ok) constant for NHWC broadcast.
  std::vector<float> inv_std_host(Ok), neg_mean_host(Ok);
  for (int o = 0; o < Ok; ++o) {
    inv_std_host[o] = 1.0f / std::sqrt(var->data[o] + kBNEpsilon);
    neg_mean_host[o] = -mean->data[o];
  }
  NSArray* bn_shape = @[@1, @1, @1, @(Ok)];
  MPSGraphTensor* mean_t = ConstFloat32(g, mean->data, bn_shape,
      [name stringByAppendingString:@"_bn_mean"]);
  MPSGraphTensor* inv_std_t = ConstFloat32(g, inv_std_host.data(), bn_shape,
      [name stringByAppendingString:@"_bn_inv_std"]);
  MPSGraphTensor* beta_t = ConstFloat32(g, beta->data, bn_shape,
      [name stringByAppendingString:@"_bn_beta"]);

  MPSGraphTensor* centered =
      [g subtractionWithPrimaryTensor:conv secondaryTensor:mean_t
                                  name:[name stringByAppendingString:@"_bn_sub"]];
  MPSGraphTensor* scaled =
      [g multiplicationWithPrimaryTensor:centered secondaryTensor:inv_std_t
                                    name:[name stringByAppendingString:@"_bn_mul"]];
  MPSGraphTensor* shifted =
      [g additionWithPrimaryTensor:scaled secondaryTensor:beta_t
                              name:[name stringByAppendingString:@"_bn_add"]];
  return [g reLUWithTensor:shifted
                       name:[name stringByAppendingString:@"_r"]];
}

// Conv-BN-ReLU: emits the fused conv + bias + relu.
MPSGraphTensor* CBR(MPSGraph* g, MPSGraphTensor* x,
                    const DvwWeights& dvw,
                    int conv_n, int bn_n,
                    int stride_y, int stride_x,
                    bool same_padding,
                    NSString* name) {
  if (g_unfold_bn_for_graph) {
    return CBRUnfolded(g, x, dvw, conv_n, bn_n, stride_y, stride_x,
                        same_padding, name);
  }
  FusedConv fc = FoldConvBn(dvw, conv_n, bn_n);
  if (fc.weights_hwio.empty()) return nullptr;
  MPSGraphTensor* y = AddConv(g, x, fc, stride_y, stride_x, same_padding, name);
  return [g reLUWithTensor:y name:[name stringByAppendingString:@"_r"]];
}

// AvgPool 3×3 + CBR.
MPSGraphTensor* AvgCBR(MPSGraph* g, MPSGraphTensor* x,
                       const DvwWeights& dvw,
                       int conv_n, int bn_n, NSString* name) {
  MPSGraphPooling2DOpDescriptor* pdesc =
      [MPSGraphPooling2DOpDescriptor
          descriptorWithKernelWidth:3
                       kernelHeight:3
                          strideInX:1
                          strideInY:1
                       paddingStyle:MPSGraphPaddingStyleTF_SAME
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC];
  // Keras AvgPool2D / DeepVariant Inception-v3 default is
  // count_include_pad=False (i.e. divide by the number of *real* kernel
  // positions, not by kernel area). MPSGraph defaults to YES, so we
  // override.
  pdesc.includeZeroPadToAverage = NO;
  MPSGraphTensor* p = [g avgPooling2DWithSourceTensor:x
                                            descriptor:pdesc
                                                  name:[name stringByAppendingString:@"_ap"]];
  return CBR(g, p, dvw, conv_n, bn_n, 1, 1, true, name);
}

MPSGraphTensor* MaxPool3x3s2Valid(MPSGraph* g, MPSGraphTensor* x,
                                  NSString* name) {
  MPSGraphPooling2DOpDescriptor* pdesc =
      [MPSGraphPooling2DOpDescriptor
          descriptorWithKernelWidth:3
                       kernelHeight:3
                          strideInX:2
                          strideInY:2
                       paddingStyle:MPSGraphPaddingStyleExplicit
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC];
  pdesc.paddingLeft = 0;
  pdesc.paddingRight = 0;
  pdesc.paddingTop = 0;
  pdesc.paddingBottom = 0;
  return [g maxPooling2DWithSourceTensor:x
                              descriptor:pdesc
                                    name:name];
}

// Inception blocks — direct ports from inception_v3_mil.py.

MPSGraphTensor* Mixed_5b(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=5 (branch1x1)
  MPSGraphTensor* b1 = CBR(g, x, d, 16, 20, 1, 1, true, @"5b_1");
  // M=6 (branch5x5 reduce)
  MPSGraphTensor* b5 = CBR(g, x, d, 12, 14, 1, 1, true, @"5b_5a");
  // M=7 (branch5x5)
  b5 = CBR(g, b5, d, 17, 21, 1, 1, true, @"5b_5b");
  // M=8 (branch3x3dbl reduce)
  MPSGraphTensor* b3 = CBR(g, x, d, 10, 11, 1, 1, true, @"5b_3a");
  // M=9 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 13, 15, 1, 1, true, @"5b_3b");
  // M=10 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 18, 22, 1, 1, true, @"5b_3c");
  // M=11 (branchpool 1×1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 19, 23, @"5b_p");
  return [g concatTensors:@[b1, b5, b3, bp] dimension:3 name:@"5b"];
}

MPSGraphTensor* Mixed_5c(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=12 (branch1x1)
  MPSGraphTensor* b1 = CBR(g, x, d, 30, 34, 1, 1, true, @"5c_1");
  // M=13 (branch5x5 reduce)
  MPSGraphTensor* b5 = CBR(g, x, d, 26, 28, 1, 1, true, @"5c_5a");
  // M=14 (branch5x5)
  b5 = CBR(g, b5, d, 31, 35, 1, 1, true, @"5c_5b");
  // M=15 (branch3x3dbl reduce)
  MPSGraphTensor* b3 = CBR(g, x, d, 24, 25, 1, 1, true, @"5c_3a");
  // M=16 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 27, 29, 1, 1, true, @"5c_3b");
  // M=17 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 32, 36, 1, 1, true, @"5c_3c");
  // M=18 (branchpool 1×1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 33, 37, @"5c_p");
  return [g concatTensors:@[b1, b5, b3, bp] dimension:3 name:@"5c"];
}

MPSGraphTensor* Mixed_5d(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=19 (branch1x1)
  MPSGraphTensor* b1 = CBR(g, x, d, 44, 48, 1, 1, true, @"5d_1");
  // M=20 (branch5x5 reduce)
  MPSGraphTensor* b5 = CBR(g, x, d, 40, 42, 1, 1, true, @"5d_5a");
  // M=21 (branch5x5)
  b5 = CBR(g, b5, d, 45, 49, 1, 1, true, @"5d_5b");
  // M=22 (branch3x3dbl reduce)
  MPSGraphTensor* b3 = CBR(g, x, d, 38, 39, 1, 1, true, @"5d_3a");
  // M=23 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 41, 43, 1, 1, true, @"5d_3b");
  // M=24 (branch3x3dbl 3×3)
  b3 = CBR(g, b3, d, 46, 50, 1, 1, true, @"5d_3c");
  // M=25 (branchpool 1×1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 47, 51, @"5d_p");
  return [g concatTensors:@[b1, b5, b3, bp] dimension:3 name:@"5d"];
}

MPSGraphTensor* Mixed_6a(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=26 (branch3x3 stride 2 valid)
  MPSGraphTensor* b3 = CBR(g, x, d, 56, 58, 2, 2, false, @"6a_3");
  // M=27 (branch3x3dbl_a 1×1)
  MPSGraphTensor* bd = CBR(g, x, d, 52, 53, 1, 1, true, @"6a_da");
  // M=28 (branch3x3dbl_b 3×3)
  bd = CBR(g, bd, d, 54, 55, 1, 1, true, @"6a_db");
  // M=29 (branch3x3dbl_c 3×3 stride 2 valid)
  bd = CBR(g, bd, d, 57, 59, 2, 2, false, @"6a_dc");
  MPSGraphTensor* bp = MaxPool3x3s2Valid(g, x, @"6a_mp");
  return [g concatTensors:@[b3, bd, bp] dimension:3 name:@"6a"];
}

MPSGraphTensor* Mixed_6b(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=30 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 72, 76, 1, 1, true, @"6b_1");
  // M=31 (b7_a, 1×1 reduce)
  MPSGraphTensor* b7a = CBR(g, x, d, 64, 66, 1, 1, true, @"6b_7aa");
  // M=32 (b7_b, 1×7)
  b7a = CBR(g, b7a, d, 68, 70, 1, 1, true, @"6b_7ab");
  // M=33 (b7_c, 7×1)
  b7a = CBR(g, b7a, d, 73, 77, 1, 1, true, @"6b_7ac");
  // M=34 (b7dbl_a, 1×1 reduce)
  MPSGraphTensor* b7b = CBR(g, x, d, 60, 61, 1, 1, true, @"6b_7ba");
  // M=35 (b7dbl_b, 7×1)
  b7b = CBR(g, b7b, d, 62, 63, 1, 1, true, @"6b_7bb");
  // M=36 (b7dbl_c, 1×7)
  b7b = CBR(g, b7b, d, 65, 67, 1, 1, true, @"6b_7bc");
  // M=37 (b7dbl_d, 7×1)
  b7b = CBR(g, b7b, d, 69, 71, 1, 1, true, @"6b_7bd");
  // M=38 (b7dbl_e, 1×7)
  b7b = CBR(g, b7b, d, 74, 78, 1, 1, true, @"6b_7be");
  // M=39 (bp_1x1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 75, 79, @"6b_p");
  return [g concatTensors:@[b1, b7a, b7b, bp] dimension:3 name:@"6b"];
}

MPSGraphTensor* Mixed_6c(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=40 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 92, 96, 1, 1, true, @"6c_1");
  // M=41 (b7_a, 1×1 reduce)
  MPSGraphTensor* b7a = CBR(g, x, d, 84, 86, 1, 1, true, @"6c_7aa");
  // M=42 (b7_b, 1×7)
  b7a = CBR(g, b7a, d, 88, 90, 1, 1, true, @"6c_7ab");
  // M=43 (b7_c, 7×1)
  b7a = CBR(g, b7a, d, 93, 97, 1, 1, true, @"6c_7ac");
  // M=44 (b7dbl_a, 1×1 reduce)
  MPSGraphTensor* b7b = CBR(g, x, d, 80, 81, 1, 1, true, @"6c_7ba");
  // M=45 (b7dbl_b, 7×1)
  b7b = CBR(g, b7b, d, 82, 83, 1, 1, true, @"6c_7bb");
  // M=46 (b7dbl_c, 1×7)
  b7b = CBR(g, b7b, d, 85, 87, 1, 1, true, @"6c_7bc");
  // M=47 (b7dbl_d, 7×1)
  b7b = CBR(g, b7b, d, 89, 91, 1, 1, true, @"6c_7bd");
  // M=48 (b7dbl_e, 1×7)
  b7b = CBR(g, b7b, d, 94, 98, 1, 1, true, @"6c_7be");
  // M=49 (bp_1x1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 95, 99, @"6c_p");
  return [g concatTensors:@[b1, b7a, b7b, bp] dimension:3 name:@"6c"];
}

MPSGraphTensor* Mixed_6d(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=50 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 112, 116, 1, 1, true, @"6d_1");
  // M=51 (b7_a, 1×1 reduce)
  MPSGraphTensor* b7a = CBR(g, x, d, 104, 106, 1, 1, true, @"6d_7aa");
  // M=52 (b7_b, 1×7)
  b7a = CBR(g, b7a, d, 108, 110, 1, 1, true, @"6d_7ab");
  // M=53 (b7_c, 7×1)
  b7a = CBR(g, b7a, d, 113, 117, 1, 1, true, @"6d_7ac");
  // M=54 (b7dbl_a, 1×1 reduce)
  MPSGraphTensor* b7b = CBR(g, x, d, 100, 101, 1, 1, true, @"6d_7ba");
  // M=55 (b7dbl_b, 7×1)
  b7b = CBR(g, b7b, d, 102, 103, 1, 1, true, @"6d_7bb");
  // M=56 (b7dbl_c, 1×7)
  b7b = CBR(g, b7b, d, 105, 107, 1, 1, true, @"6d_7bc");
  // M=57 (b7dbl_d, 7×1)
  b7b = CBR(g, b7b, d, 109, 111, 1, 1, true, @"6d_7bd");
  // M=58 (b7dbl_e, 1×7)
  b7b = CBR(g, b7b, d, 114, 118, 1, 1, true, @"6d_7be");
  // M=59 (bp_1x1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 115, 119, @"6d_p");
  return [g concatTensors:@[b1, b7a, b7b, bp] dimension:3 name:@"6d"];
}

MPSGraphTensor* Mixed_6e(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=60 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 132, 136, 1, 1, true, @"6e_1");
  // M=61 (b7_a, 1×1 reduce)
  MPSGraphTensor* b7a = CBR(g, x, d, 124, 126, 1, 1, true, @"6e_7aa");
  // M=62 (b7_b, 1×7)
  b7a = CBR(g, b7a, d, 128, 130, 1, 1, true, @"6e_7ab");
  // M=63 (b7_c, 7×1)
  b7a = CBR(g, b7a, d, 133, 137, 1, 1, true, @"6e_7ac");
  // M=64 (b7dbl_a, 1×1 reduce)
  MPSGraphTensor* b7b = CBR(g, x, d, 120, 121, 1, 1, true, @"6e_7ba");
  // M=65 (b7dbl_b, 7×1)
  b7b = CBR(g, b7b, d, 122, 123, 1, 1, true, @"6e_7bb");
  // M=66 (b7dbl_c, 1×7)
  b7b = CBR(g, b7b, d, 125, 127, 1, 1, true, @"6e_7bc");
  // M=67 (b7dbl_d, 7×1)
  b7b = CBR(g, b7b, d, 129, 131, 1, 1, true, @"6e_7bd");
  // M=68 (b7dbl_e, 1×7)
  b7b = CBR(g, b7b, d, 134, 138, 1, 1, true, @"6e_7be");
  // M=69 (bp_1x1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 135, 139, @"6e_p");
  return [g concatTensors:@[b1, b7a, b7b, bp] dimension:3 name:@"6e"];
}

MPSGraphTensor* Mixed_7a(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=70 (b3_a, 1×1)
  MPSGraphTensor* b3 = CBR(g, x, d, 144, 146, 1, 1, true, @"7a_3a");
  // M=71 (b3_b, 3×3 stride 2 valid)
  b3 = CBR(g, b3, d, 148, 150, 2, 2, false, @"7a_3b");
  // M=72 (b7_a, 1×1)
  MPSGraphTensor* b7 = CBR(g, x, d, 140, 141, 1, 1, true, @"7a_7a");
  // M=73 (b7_b, 1×7)
  b7 = CBR(g, b7, d, 142, 143, 1, 1, true, @"7a_7b");
  // M=74 (b7_c, 7×1)
  b7 = CBR(g, b7, d, 145, 147, 1, 1, true, @"7a_7c");
  // M=75 (b7_d, 3×3 stride 2 valid)
  b7 = CBR(g, b7, d, 149, 151, 2, 2, false, @"7a_7d");
  MPSGraphTensor* bp = MaxPool3x3s2Valid(g, x, @"7a_mp");
  return [g concatTensors:@[b3, b7, bp] dimension:3 name:@"7a"];
}

MPSGraphTensor* Mixed_7b(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=76 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 162, 168, 1, 1, true, @"7b_1");
  // M=77 (b3 reduce 1×1)
  MPSGraphTensor* b3a = CBR(g, x, d, 154, 156, 1, 1, true, @"7b_3aa");
  // M=78 (b3 1×3)
  MPSGraphTensor* b3a_1x3 = CBR(g, b3a, d, 158, 163, 1, 1, true, @"7b_3a1x3");
  // M=79 (b3 3×1)
  MPSGraphTensor* b3a_3x1 = CBR(g, b3a, d, 159, 164, 1, 1, true, @"7b_3a3x1");
  // M=80 (b3dbl reduce 1×1)
  MPSGraphTensor* b3b = CBR(g, x, d, 152, 153, 1, 1, true, @"7b_3ba");
  // M=81 (b3dbl 3×3)
  b3b = CBR(g, b3b, d, 155, 157, 1, 1, true, @"7b_3bb");
  // M=82 (b3dbl 1×3)
  MPSGraphTensor* b3b_1x3 = CBR(g, b3b, d, 160, 165, 1, 1, true, @"7b_3b1x3");
  // M=83 (b3dbl 3×1)
  MPSGraphTensor* b3b_3x1 = CBR(g, b3b, d, 161, 166, 1, 1, true, @"7b_3b3x1");
  // M=84 (bp 1×1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 167, 169, @"7b_p");
  return [g concatTensors:@[b1, b3a_1x3, b3a_3x1, b3b_1x3, b3b_3x1, bp] dimension:3 name:@"7b"];
}

MPSGraphTensor* Mixed_7c(MPSGraph* g, MPSGraphTensor* x,
                         const DvwWeights& d) {
  // M=85 (b1, 1×1)
  MPSGraphTensor* b1 = CBR(g, x, d, 180, 186, 1, 1, true, @"7c_1");
  // M=86 (b3 reduce 1×1)
  MPSGraphTensor* b3a = CBR(g, x, d, 172, 174, 1, 1, true, @"7c_3aa");
  // M=87 (b3 1×3)
  MPSGraphTensor* b3a_1x3 = CBR(g, b3a, d, 176, 181, 1, 1, true, @"7c_3a1x3");
  // M=88 (b3 3×1)
  MPSGraphTensor* b3a_3x1 = CBR(g, b3a, d, 177, 182, 1, 1, true, @"7c_3a3x1");
  // M=89 (b3dbl reduce 1×1)
  MPSGraphTensor* b3b = CBR(g, x, d, 170, 171, 1, 1, true, @"7c_3ba");
  // M=90 (b3dbl 3×3)
  b3b = CBR(g, b3b, d, 173, 175, 1, 1, true, @"7c_3bb");
  // M=91 (b3dbl 1×3)
  MPSGraphTensor* b3b_1x3 = CBR(g, b3b, d, 178, 183, 1, 1, true, @"7c_3b1x3");
  // M=92 (b3dbl 3×1)
  MPSGraphTensor* b3b_3x1 = CBR(g, b3b, d, 179, 184, 1, 1, true, @"7c_3b3x1");
  // M=93 (bp 1×1)
  MPSGraphTensor* bp = AvgCBR(g, x, d, 185, 187, @"7c_p");
  return [g concatTensors:@[b1, b3a_1x3, b3a_3x1, b3b_1x3, b3b_3x1, bp] dimension:3 name:@"7c"];
}


}  // namespace

// ---------------------------------------------------------------------------
// Impl: holds device, queue, graph, and the cached executable.
// ---------------------------------------------------------------------------

// Deterministic-stage params. Either a Conv2D (with weights+bias) or a
// MaxPool2D (no weights). Populated when DV_METAL_DET_LAYERS names a
// stage. Each stage takes the previous stage's output as input and
// writes its own output to a dst buffer at Predict time.
struct DetLayer {
  std::string tap_name;     // "stem_s1a" / "stem_mp3a" — output tap name
  enum Kind { kConv, kMaxPool } kind = kConv;
  // Common geometry:
  int C_in = 0, C_out = 0;
  int H_in = 0, W_in = 0;
  int H_out = 0, W_out = 0;
  // Conv-specific (kind == kConv):
  ConvDesc conv_desc{};
  id<MTLBuffer> weights_buf = nil;
  id<MTLBuffer> bias_buf = nil;
  // Phase 5.5f — unfolded conv→BN→ReLU. When `use_unfolded_bn` is set,
  // weights_buf holds RAW kernel HWIO (no inv_std scaling) and bias_buf
  // is an all-zero buffer. The conv is encoded with relu=false; a
  // separate MetalBnRelu pass consumes its output using the BN params
  // below and applies ReLU. Bit-match measurement (Phase 5.5f Day 1)
  // shows this path matches TF/oneDNN to ±2 ULP per element vs ±93 ULP
  // for the folded path.
  bool use_unfolded_bn = false;
  id<MTLBuffer> bn_mean_buf = nil;
  id<MTLBuffer> bn_var_buf = nil;
  id<MTLBuffer> bn_beta_buf = nil;
  id<MTLBuffer> bn_inter_buf = nil;     // post-conv pre-BN intermediate
  // MaxPool-specific (kind == kMaxPool):
  MaxPoolDesc pool_desc{};
  // The MPSGraph "post-graph" — only populated on the LAST det stage in
  // the chain. Takes that stage's output as placeholder, runs through
  // gap.
  MPSGraph* post_graph = nil;
  MPSGraphTensor* post_input = nil;
  MPSGraphTensor* post_output = nil;
};

struct MetalInception::Impl {
  std::unique_ptr<DvwWeights> weights;
  // Input shape parameters: NHWC.
  // WGS: H=100 W=221 C=7. DeepTrio WGS: H=140 W=221 C=7.
  // PacBio: H=100 W=147 C=10. ONT: H=100 W=199 C=10.
  // Somatic PacBio TN: H=200 W=147 C=9.
  int input_height = 100;
  int input_width = 221;
  int input_channels = 7;
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  MPSGraph* graph = nil;
  MPSGraphTensor* input = nil;
  MPSGraphTensor* output = nil;
  // ── DV_METAL_GPU_FINALIZE=1 (default off) ─────────────────────────────
  // When set, append the (2048→3) dense + softmax to the graph so
  // Predict() returns probabilities (B,3) instead of features (B,2048).
  // Bypasses the BnnsFinalize CPU step. Outputs are GPU softmax via
  // MPSGraph's parallel reduction (different rounding from BNNS-CPU
  // sequential), so a per-chip drift on the order of ~1 ULP at the
  // softmax may differ from the BNNS-CPU baseline.
  bool gpu_finalize = false;
  int output_dim = 2048;  // 2048 by default; 3 with gpu_finalize=true
  // Named taps for debugging — keyed by stage name. Populated as the
  // graph is built, so PredictAtTap() can request a specific stage's
  // output.
  NSMutableDictionary<NSString*, MPSGraphTensor*>* taps = nil;
  // Compiled executable with optimizationLevel=Level0 (GPU-only,
  // no ANE placement pass). Lazily filled per tap on first request.
  MPSGraphCompilationDescriptor* compileDesc = nil;
  NSMutableDictionary<NSString*, MPSGraphExecutable*>* execCache = nil;
  // Ordered list of tap tensors compiled into the gap executable
  // (used for the full-network forward — Predict()).
  int feature_dim = 2048;

  // ── Phase 5.5c — deterministic-reduction-order stem path ──────────
  // When det_layers is empty, Predict() takes the original full-graph
  // path. When non-empty, Predict() runs: det kernels in chain → post-
  // graph → gap. Det stages are a contiguous prefix of the network's
  // first 7 layers (s1a, s2a, s2b, mp3a, s3b, s4a, mp5a).
  std::unique_ptr<MetalConvSerial> conv_serial;
  std::unique_ptr<MetalMaxPool> max_pool;
  std::unique_ptr<MetalBnRelu> bn_relu;     // Phase 5.5f, lazy-init
  std::vector<DetLayer> det_layers;
  // Phase 8 / Tier 6.0 — full-network det Inception path. When non-empty,
  // Predict() runs det stem chain → det_blocks chain → global_avg_pool →
  // output, completely bypassing MPSGraph on the conv path. Bit-deterministic
  // across runs/chips. Activated by DV_METAL_SERIAL_FULL=1 env var.
  std::vector<DetMixedBlock> det_blocks;
  std::unique_ptr<MetalAvgPool> avg_pool;
  std::unique_ptr<MetalConcat> concat;
  std::unique_ptr<MetalGlobalAvgPool> gap_pool;
  id<MTLBuffer> gap_out_buf = nil;          // (max_B, 2048) FP32 — output of global avg pool
  // Cached post-graph executable per batch size (the post-graph itself is
  // already in det_layers[0].post_graph).
  NSMutableDictionary<NSNumber*, MPSGraphExecutable*>* post_exec_cache = nil;
};

MetalInception::MetalInception() : impl_(std::make_unique<Impl>()) {}
MetalInception::~MetalInception() = default;

int MetalInception::FeatureDim() const {
  // Returns the per-example output dimension Predict() writes:
  //   - default path: feature_dim (2048 for standard Inception-v3, det
  //     chain may override via Mixed_7c output channel count)
  //   - DV_METAL_GPU_FINALIZE=1: 3 (post-softmax probabilities)
  if (!impl_) return 0;
  return impl_->gpu_finalize ? impl_->output_dim : impl_->feature_dim;
}

bool MetalInception::IsGpuFinalize() const {
  return impl_ && impl_->gpu_finalize;
}

std::unique_ptr<MetalInception> MetalInception::Create(
    const std::string& dvw_path,
    int input_height,
    int input_channels,
    int input_width) {
  auto self = std::unique_ptr<MetalInception>(new MetalInception());
  auto& I = *self->impl_;
  I.input_height = input_height;
  I.input_width = input_width;
  I.input_channels = input_channels;

  // Phase 5.5f: read DV_METAL_UNFOLDED_BN before any graph stages are
  // built so CBR()/CBRUnfolded() dispatch correctly throughout. This
  // global flag is also reused by the det-path env-var check below.
  {
    const char* env = std::getenv("DV_METAL_UNFOLDED_BN");
    g_unfold_bn_for_graph =
        (env && std::string(env) != "0" && std::string(env) != "false");
    if (g_unfold_bn_for_graph) {
      LOG(INFO) << "Phase 5.5f: unfolded conv→BN→ReLU active for full graph "
                << "(every CBR call uses raw conv + primitive BN ops)";
    }
  }

  I.weights = DvwWeights::Open(dvw_path);
  if (!I.weights) {
    LOG(ERROR) << "MetalInception::Create: cannot open " << dvw_path;
    return nullptr;
  }

  I.device = MTLCreateSystemDefaultDevice();
  if (!I.device) {
    LOG(ERROR) << "MetalInception::Create: no Metal device available";
    return nullptr;
  }
  I.queue = [I.device newCommandQueue];
  if (!I.queue) {
    LOG(ERROR) << "MetalInception::Create: failed to create command queue";
    return nullptr;
  }

  I.graph = [MPSGraph new];
  // Compilation descriptor: optimizationLevel=Level0 disables the
  // "placement pass dispatching across NeuralEngine and CPU along
  // with the GPU" (per MPSGraph.h). Default Level1 silently picks
  // mixed-precision paths (e.g. FP16 Winograd intermediates) and
  // off-GPU placements for ops where it thinks it's safe — which
  // produces channel-permuted output for our FP32 Inception-v3 conv.
  // Level0 forces GPU-only, full-precision execution at the cost of
  // some perf optimisations.
  // DV_METAL_FP16=1 opts into reduced-precision inference: MPSGraph is
  // allowed to use FP16 Winograd intermediates + operand conversion and the
  // Level1 perf optimisations. This is faster on Apple Silicon but is NOT
  // bit-identical to the FP32 path and changes genotype likelihoods at the
  // last digits — gate any default flip behind a GIAB concordance check.
  // Default (unset) keeps the contractual full-FP32 behaviour unchanged.
  const char* fp16_env = std::getenv("DV_METAL_FP16");
  const bool fp16 = fp16_env && fp16_env[0] == '1';
  I.compileDesc = [MPSGraphCompilationDescriptor new];
  I.compileDesc.optimizationLevel =
      fp16 ? MPSGraphOptimizationLevel1 : MPSGraphOptimizationLevel0;
  I.compileDesc.waitForCompilationCompletion = YES;
  // macOS 26+: control reduced-precision fast-math paths (FP16 Winograd
  // intermediates, TF32 (19-bit) operand conversion). Default is `None` already
  // — set explicitly to make the behaviour contractually visible and logged.
  if (@available(macOS 26.0, iOS 26.0, *)) {
    I.compileDesc.reducedPrecisionFastMath =
        fp16 ? MPSGraphReducedPrecisionFastMathDefault
             : MPSGraphReducedPrecisionFastMathNone;
    LOG(INFO) << "MPSGraph: reducedPrecisionFastMath="
              << (fp16 ? "Default (FP16 fast-math allowed)"
                       : "None (full FP32)");
  }
  if (fp16) {
    LOG(INFO) << "MetalInception: DV_METAL_FP16=1 — reduced-precision "
                 "inference (Level1, FP16 fast-math). NOT bit-identical to "
                 "FP32; validate accuracy before relying on it.";
  }
  I.execCache = [NSMutableDictionary dictionary];
  I.taps = [NSMutableDictionary dictionary];
  // Variable batch dimension. -1 means "any" in MPSGraph shape spec.
  // Height/channels parameterized at construction (WGS=100/7, trio=140/7).
  I.input = [I.graph placeholderWithShape:@[@-1,
                                              @(I.input_height),
                                              @(I.input_width),
                                              @(I.input_channels)]
                                  dataType:MPSDataTypeFloat32
                                      name:@"input_nhwc"];
  // Stay in NHWC throughout — TF native layout. (Earlier OIHW/NCHW path
  // produced channel-permuted output despite a hand-rolled transpose
  // matching TF; switching to NHWC end-to-end resolved it.)
  MPSGraphTensor* x = I.input;
  I.taps[@"input_nchw"] = x;  // tap kept under the old name; layout = NHWC now

  // Stem
  x = CBR(I.graph, x, *I.weights, 0, 1, 2, 2, false, @"s1a");
  if (!x) return nullptr;
  I.taps[@"stem_s1a"] = x;
  x = CBR(I.graph, x, *I.weights, 2, 3, 1, 1, false, @"s2a");
  I.taps[@"stem_s2a"] = x;
  x = CBR(I.graph, x, *I.weights, 4, 5, 1, 1, true,  @"s2b");
  I.taps[@"stem_s2b"] = x;
  x = MaxPool3x3s2Valid(I.graph, x, @"mp3a");
  I.taps[@"stem_mp3a"] = x;
  x = CBR(I.graph, x, *I.weights, 6, 7, 1, 1, false, @"s3b");
  I.taps[@"stem_s3b"] = x;
  x = CBR(I.graph, x, *I.weights, 8, 9, 1, 1, false, @"s4a");
  I.taps[@"stem_s4a"] = x;
  x = MaxPool3x3s2Valid(I.graph, x, @"mp5a");
  I.taps[@"stem_mp5a"] = x;

  // InceptionA
  x = Mixed_5b(I.graph, x, *I.weights); I.taps[@"5b"] = x;
  x = Mixed_5c(I.graph, x, *I.weights); I.taps[@"5c"] = x;
  x = Mixed_5d(I.graph, x, *I.weights); I.taps[@"5d"] = x;
  // Reduction-A
  x = Mixed_6a(I.graph, x, *I.weights); I.taps[@"6a"] = x;
  // InceptionB
  x = Mixed_6b(I.graph, x, *I.weights); I.taps[@"6b"] = x;
  x = Mixed_6c(I.graph, x, *I.weights); I.taps[@"6c"] = x;
  x = Mixed_6d(I.graph, x, *I.weights); I.taps[@"6d"] = x;
  x = Mixed_6e(I.graph, x, *I.weights); I.taps[@"6e"] = x;
  // Reduction-B
  x = Mixed_7a(I.graph, x, *I.weights); I.taps[@"7a"] = x;
  // InceptionC
  x = Mixed_7b(I.graph, x, *I.weights); I.taps[@"7b"] = x;
  x = Mixed_7c(I.graph, x, *I.weights); I.taps[@"7c"] = x;

  // Global avg pool over (H, W) → (N, 2048, 1, 1)
  x = [I.graph meanOfTensor:x axes:@[@1, @2] name:@"gap"];
  // Reshape to (N, 2048)
  x = [I.graph reshapeTensor:x withShape:@[@-1, @2048] name:@"squeeze"];
  I.taps[@"gap"] = x;

  I.output = x;
  I.output_dim = 2048;

  // ── DV_METAL_GPU_FINALIZE=1: append dense (2048→3) + softmax ──────────
  // The terminal classifier head moves from BnnsFinalize (sequential
  // FP32 on CPU) to MPSGraph (parallel reduction on GPU). One less
  // host-device sync per batch; functional equivalence on chr20 to be
  // verified by FILTER-class diff vs the BNNS-CPU baseline.
  {
    const char* gf_env = std::getenv("DV_METAL_GPU_FINALIZE");
    if (gf_env && *gf_env && std::string(gf_env) != "0") {
      const auto* k = I.weights->Get(
          "layer_with_weights-188/kernel/.ATTRIBUTES/VARIABLE_VALUE");
      const auto* b = I.weights->Get(
          "layer_with_weights-188/bias/.ATTRIBUTES/VARIABLE_VALUE");
      if (!k || !b || k->shape.size() != 2u || b->shape.size() != 1u ||
          k->shape[0] != 2048u || k->shape[1] != 3u || b->shape[0] != 3u) {
        LOG(ERROR) << "MetalInception::Create: DV_METAL_GPU_FINALIZE=1 set "
                      "but layer-188 weights are missing or wrong shape "
                   << "(expected kernel (2048,3) + bias (3,))";
        return nullptr;
      }
      MPSGraphTensor* W = ConstFloat32(I.graph, k->data,
                                       @[@2048, @3], @"finalize_w");
      MPSGraphTensor* B = ConstFloat32(I.graph, b->data,
                                       @[@3], @"finalize_b");
      // logits = features (N,2048) · W (2048,3) → (N,3)
      MPSGraphTensor* logits =
          [I.graph matrixMultiplicationWithPrimaryTensor:x
                                          secondaryTensor:W
                                                    name:@"finalize_matmul"];
      // bias broadcast (3,) → (1,3)
      MPSGraphTensor* B_r = [I.graph reshapeTensor:B
                                          withShape:@[@1, @3]
                                               name:@"finalize_b_r"];
      logits = [I.graph additionWithPrimaryTensor:logits
                                  secondaryTensor:B_r
                                             name:@"finalize_logits"];
      // softmax along channel axis (axis 1 in (N,3))
      MPSGraphTensor* probs = [I.graph softMaxWithTensor:logits
                                                    axis:1
                                                    name:@"finalize_softmax"];
      I.taps[@"probs"] = probs;
      I.output = probs;
      I.output_dim = 3;
      I.gpu_finalize = true;
      LOG(INFO) << "MetalInception: DV_METAL_GPU_FINALIZE=1 — "
                << "graph outputs (N,3) probs (BnnsFinalize bypassed)";
    }
  }

  // ── Phase 5.5c: optional deterministic kernel for stem_s1a ─────────
  // DV_METAL_DET_LAYERS=stem_s1a triggers a parallel inference path
  // where the first conv (CBR(0,1) stride 2-2 valid 7→32) is replaced
  // by a deterministic-reduction-order Metal compute kernel and the
  // network from stem_s2a through gap is run via a separate MPSGraph
  // that takes the s1a output as a placeholder. Other layer names are
  // ignored for now (extension to additional CBR layers is a follow-up).
  const char* det_env = std::getenv("DV_METAL_DET_LAYERS");
  std::set<std::string> det_set;
  if (det_env && *det_env) {
    for (absl::string_view s : absl::StrSplit(det_env, ',')) {
      if (!s.empty()) det_set.emplace(s.data(), s.size());
    }
  }
  // Det layers must form a contiguous chain starting at the head of
  // the stem. Supported stages (in order):
  //   stem_s1a (CBR 0,1) stride 2 VALID 7→32  : (B,100,221,7) → (B,49,110,32)
  //   stem_s2a (CBR 2,3) stride 1 VALID 32→32 : (B,49,110,32) → (B,47,108,32)
  //   stem_s2b (CBR 4,5) stride 1 SAME  32→64 : (B,47,108,32) → (B,47,108,64)
  //   stem_mp3a maxpool 3×3 stride 2 VALID    : (B,47,108,64) → (B,23,53,64)
  //   stem_s3b (CBR 6,7) stride 1 VALID 64→80 : (B,23,53,64) → (B,23,53,80)
  //   stem_s4a (CBR 8,9) stride 1 VALID 80→192: (B,23,53,80) → (B,21,51,192)
  //   stem_mp5a maxpool 3×3 stride 2 VALID    : (B,21,51,192)→ (B,10,25,192)
  // Convenience: DV_METAL_DET_LAYERS=stem expands to all 7.
  enum StemKind { kSCBR, kSPool };
  struct StemStage { const char* tap; StemKind kind;
                      // CBR-only:
                      int conv; int bn; int sy; int sx; bool same;
                      // Geometry (always set):
                      int C_in; int C_out;
                      int H_in; int W_in; int H_out; int W_out; };
  static const StemStage kStemChain[] = {
      {"stem_s1a",  kSCBR, 0, 1, 2, 2, false,    7,  32, 100, 221, 49, 110},
      {"stem_s2a",  kSCBR, 2, 3, 1, 1, false,   32,  32,  49, 110, 47, 108},
      {"stem_s2b",  kSCBR, 4, 5, 1, 1, true,    32,  64,  47, 108, 47, 108},
      {"stem_mp3a", kSPool, 0, 0, 2, 2, false,  64,  64,  47, 108, 23,  53},
      {"stem_s3b",  kSCBR, 6, 7, 1, 1, false,   64,  80,  23,  53, 23,  53},
      {"stem_s4a",  kSCBR, 8, 9, 1, 1, false,   80, 192,  23,  53, 21,  51},
      {"stem_mp5a", kSPool, 0, 0, 2, 2, false, 192, 192,  21,  51, 10,  25},
  };
  // "stem" alias enables ALL 7 stages.
  if (det_set.count("stem") > 0) {
    for (const auto& s : kStemChain) det_set.insert(s.tap);
  }
  int chain_len = 0;
  for (const auto& s : kStemChain) {
    if (det_set.count(s.tap) > 0) {
      if (chain_len == &s - kStemChain) ++chain_len;
      else {
        LOG(ERROR)
            << "DV_METAL_DET_LAYERS: must be contiguous from stem_s1a; "
               "got non-contiguous set";
        return nullptr;
      }
    }
  }
  if (chain_len > 0) {
    LOG(INFO) << "Metal det path: " << chain_len
              << " stem stage(s) → deterministic kernel chain";

    I.conv_serial = MetalConvSerial::Create();
    I.max_pool = MetalMaxPool::Create();
    if (!I.conv_serial || !I.max_pool) {
      LOG(ERROR) << "MetalInception::Create: kernel pipeline creation failed";
      return nullptr;
    }

    // Phase 5.5f — DV_METAL_UNFOLDED_BN=1 enables the conv→BN→ReLU
    // separation that bit-matches TF/oneDNN to ±2 ULP. Without it, the
    // det path uses FoldConvBn which drifts up to 93 ULP per element.
    const char* unfolded_env = std::getenv("DV_METAL_UNFOLDED_BN");
    const bool unfolded_bn =
        (unfolded_env && std::string(unfolded_env) != "0" &&
         std::string(unfolded_env) != "false");
    if (unfolded_bn) {
      I.bn_relu = MetalBnRelu::Create();
      if (!I.bn_relu) {
        LOG(ERROR) << "MetalInception::Create: BN+ReLU pipeline failed";
        return nullptr;
      }
      LOG(INFO) << "Phase 5.5f: unfolded conv→BN→ReLU active for det path";
    }
    I.post_exec_cache = [NSMutableDictionary dictionary];

    // Build a det entry for each stage in the chain.
    for (int li = 0; li < chain_len; ++li) {
      const StemStage& s = kStemChain[li];
      DetLayer det{};
      det.tap_name = s.tap;
      det.kind = (s.kind == kSCBR) ? DetLayer::kConv : DetLayer::kMaxPool;
      det.C_in = s.C_in;
      det.C_out = s.C_out;
      det.H_in = s.H_in;
      det.W_in = s.W_in;
      det.H_out = s.H_out;
      det.W_out = s.W_out;
      if (s.kind == kSCBR) {
        if (unfolded_bn) {
          // Phase 5.5f: load raw conv kernel + raw BN params; defer fold.
          const auto* k = I.weights->Get(AttrCpp(s.conv, "kernel"));
          const auto* beta = I.weights->Get(AttrCpp(s.bn, "beta"));
          const auto* mean = I.weights->Get(AttrCpp(s.bn, "moving_mean"));
          const auto* var = I.weights->Get(AttrCpp(s.bn, "moving_variance"));
          if (!k || !beta || !mean || !var || k->shape.size() != 4u) {
            LOG(ERROR) << "MetalInception::Create: missing raw weight for "
                       << s.tap;
            return nullptr;
          }
          const int Hk = k->shape[0], Wk = k->shape[1];
          const int Ik = k->shape[2], Ok = k->shape[3];
          det.conv_desc.C_in = Ik;
          det.conv_desc.C_out = Ok;
          det.conv_desc.Kh = Hk;
          det.conv_desc.Kw = Wk;
          det.conv_desc.stride_h = s.sy;
          det.conv_desc.stride_w = s.sx;
          det.conv_desc.pad_h = s.same ? (Hk - 1) / 2 : 0;
          det.conv_desc.pad_w = s.same ? (Wk - 1) / 2 : 0;
          det.conv_desc.relu = false;   // ReLU happens after BN
          det.use_unfolded_bn = true;
          det.weights_buf =
              [I.device newBufferWithBytes:k->data
                                    length:k->n_bytes
                                   options:MTLResourceStorageModeShared];
          // All-zeros bias so conv output stays raw.
          std::vector<float> zero_bias(Ok, 0.0f);
          det.bias_buf =
              [I.device newBufferWithBytes:zero_bias.data()
                                    length:Ok * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          det.bn_mean_buf =
              [I.device newBufferWithBytes:mean->data
                                    length:Ok * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          det.bn_var_buf =
              [I.device newBufferWithBytes:var->data
                                    length:Ok * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          det.bn_beta_buf =
              [I.device newBufferWithBytes:beta->data
                                    length:Ok * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          if (!det.weights_buf || !det.bias_buf || !det.bn_mean_buf ||
              !det.bn_var_buf || !det.bn_beta_buf) {
            LOG(ERROR) << "MetalInception::Create: alloc failed for " << s.tap;
            return nullptr;
          }
        } else {
          FusedConv fc = FoldConvBn(*I.weights, s.conv, s.bn);
          if (fc.weights_hwio.empty()) {
            LOG(ERROR) << "MetalInception::Create: failed to fold "
                       << s.tap << " weights";
            return nullptr;
          }
          det.conv_desc.C_in = fc.I;
          det.conv_desc.C_out = fc.O;
          det.conv_desc.Kh = fc.H;
          det.conv_desc.Kw = fc.W;
          det.conv_desc.stride_h = s.sy;
          det.conv_desc.stride_w = s.sx;
          det.conv_desc.pad_h = s.same ? (fc.H - 1) / 2 : 0;
          det.conv_desc.pad_w = s.same ? (fc.W - 1) / 2 : 0;
          det.conv_desc.relu = true;
          det.weights_buf =
              [I.device newBufferWithBytes:fc.weights_hwio.data()
                                    length:fc.weights_hwio.size() * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          det.bias_buf =
              [I.device newBufferWithBytes:fc.bias.data()
                                    length:fc.bias.size() * sizeof(float)
                                   options:MTLResourceStorageModeShared];
          if (!det.weights_buf || !det.bias_buf) {
            LOG(ERROR) << "MetalInception::Create: alloc failed for " << s.tap;
            return nullptr;
          }
        }
      } else {
        // MaxPool: 3×3 stride-2 VALID, no learned params.
        det.pool_desc.C = s.C_in;
        det.pool_desc.Kh = 3;
        det.pool_desc.Kw = 3;
        det.pool_desc.stride_h = s.sy;
        det.pool_desc.stride_w = s.sx;
        det.pool_desc.pad_h = 0;
        det.pool_desc.pad_w = 0;
      }
      I.det_layers.push_back(std::move(det));
    }

    // Phase 8 / Tier 6.0 — full-network det path: when DV_METAL_SERIAL_FULL=1
    // is set AND the stem chain is full (s1a→mp5a) AND unfolded BN is on,
    // build all 11 Inception blocks (Mixed_5b…7c) + global avg pool. Predict()
    // will route through them, bypassing MPSGraph entirely on the conv path.
    const char* serial_full_env = std::getenv("DV_METAL_SERIAL_FULL");
    const bool serial_full =
        (serial_full_env && std::string(serial_full_env) != "0" &&
         std::string(serial_full_env) != "false");
    if (serial_full && chain_len == 7) {
      LOG(INFO) << "Phase 8/Tier 6.0: building full-network det path "
                << "(11 Inception blocks + global avg pool)";
      // Geometry input to Mixed_5b = output of stem_mp5a.
      const DetLayer& last_stem = I.det_layers.back();
      int blk_H = last_stem.H_out;
      int blk_W = last_stem.W_out;
      int blk_C = last_stem.C_out;
      // Use a generous max_B; we don't know batch size at Create time.
      // Allocations scale ~ 250 MB total across 11 blocks at B=128, fine on
      // M4 Max unified memory.
      const int max_B = 2048;
      I.avg_pool = MetalAvgPool::Create();
      I.concat = MetalConcat::Create();
      I.gap_pool = MetalGlobalAvgPool::Create();
      if (!I.avg_pool || !I.concat || !I.gap_pool) {
        LOG(ERROR) << "Tier 6.0: avg_pool/concat/gap_pool create failed";
        return nullptr;
      }
      using BuilderFn = bool(*)(id<MTLDevice>, const DvwWeights&, int, int, int, int, DetMixedBlock*);
      static const std::pair<BuilderFn, const char*> kBlockSpecs[] = {
          {BuildDetMixed5b, "5b"}, {BuildDetMixed5c, "5c"},
          {BuildDetMixed5d, "5d"}, {BuildDetMixed6a, "6a"},
          {BuildDetMixed6b, "6b"}, {BuildDetMixed6c, "6c"},
          {BuildDetMixed6d, "6d"}, {BuildDetMixed6e, "6e"},
          {BuildDetMixed7a, "7a"}, {BuildDetMixed7b, "7b"},
          {BuildDetMixed7c, "7c"},
      };
      I.det_blocks.resize(sizeof(kBlockSpecs) / sizeof(kBlockSpecs[0]));
      for (size_t i = 0; i < I.det_blocks.size(); ++i) {
        if (!kBlockSpecs[i].first(I.device, *I.weights, max_B, blk_H, blk_W, blk_C,
                                   &I.det_blocks[i])) {
          LOG(ERROR) << "Tier 6.0: build " << kBlockSpecs[i].second << " failed";
          return nullptr;
        }
        blk_H = I.det_blocks[i].H_out;
        blk_W = I.det_blocks[i].W_out;
        blk_C = I.det_blocks[i].C_out;
      }
      // Allocate gap output buffer (max_B, C_out_7c=2048).
      const size_t gap_bytes = (size_t)max_B * blk_C * sizeof(float);
      I.gap_out_buf =
          [I.device newBufferWithLength:gap_bytes
                                options:MTLResourceStorageModeShared];
      if (!I.gap_out_buf) {
        LOG(ERROR) << "Tier 6.0: gap_out_buf alloc failed";
        return nullptr;
      }
      I.feature_dim = blk_C;
      LOG(INFO) << "Phase 8/Tier 6.0: " << I.det_blocks.size()
                << " Inception blocks built; gap output dim = " << blk_C;
    }

    // Build ONE post-graph that starts after the last det stage.
    const DetLayer& last = I.det_layers.back();
    const std::string& last_tap = last.tap_name;
    DetLayer& last_mut = I.det_layers.back();
    last_mut.post_graph = [MPSGraph new];
    last_mut.post_input = [last_mut.post_graph
        placeholderWithShape:@[@-1, @(last.H_out), @(last.W_out),
                                @(last.C_out)]
                    dataType:MPSDataTypeFloat32
                        name:@"det_chain_out"];
    MPSGraphTensor* px = last_mut.post_input;
    // Append remaining stem stages (with index > last_tap's index).
    // Order: 0:s1a, 1:s2a, 2:s2b, 3:mp3a, 4:s3b, 5:s4a, 6:mp5a.
    int last_idx = -1;
    if (last_tap == "stem_s1a")  last_idx = 0;
    else if (last_tap == "stem_s2a")  last_idx = 1;
    else if (last_tap == "stem_s2b")  last_idx = 2;
    else if (last_tap == "stem_mp3a") last_idx = 3;
    else if (last_tap == "stem_s3b")  last_idx = 4;
    else if (last_tap == "stem_s4a")  last_idx = 5;
    else if (last_tap == "stem_mp5a") last_idx = 6;
    if (last_idx < 1) px = CBR(last_mut.post_graph, px, *I.weights, 2, 3, 1, 1, false, @"p_s2a");
    if (last_idx < 2) px = CBR(last_mut.post_graph, px, *I.weights, 4, 5, 1, 1, true,  @"p_s2b");
    if (last_idx < 3) px = MaxPool3x3s2Valid(last_mut.post_graph, px, @"p_mp3a");
    if (last_idx < 4) px = CBR(last_mut.post_graph, px, *I.weights, 6, 7, 1, 1, false, @"p_s3b");
    if (last_idx < 5) px = CBR(last_mut.post_graph, px, *I.weights, 8, 9, 1, 1, false, @"p_s4a");
    if (last_idx < 6) px = MaxPool3x3s2Valid(last_mut.post_graph, px, @"p_mp5a");
    px = Mixed_5b(last_mut.post_graph, px, *I.weights);
    px = Mixed_5c(last_mut.post_graph, px, *I.weights);
    px = Mixed_5d(last_mut.post_graph, px, *I.weights);
    px = Mixed_6a(last_mut.post_graph, px, *I.weights);
    px = Mixed_6b(last_mut.post_graph, px, *I.weights);
    px = Mixed_6c(last_mut.post_graph, px, *I.weights);
    px = Mixed_6d(last_mut.post_graph, px, *I.weights);
    px = Mixed_6e(last_mut.post_graph, px, *I.weights);
    px = Mixed_7a(last_mut.post_graph, px, *I.weights);
    px = Mixed_7b(last_mut.post_graph, px, *I.weights);
    px = Mixed_7c(last_mut.post_graph, px, *I.weights);
    px = [last_mut.post_graph meanOfTensor:px axes:@[@1, @2] name:@"p_gap"];
    px = [last_mut.post_graph reshapeTensor:px withShape:@[@-1, @2048]
                                       name:@"p_squeeze"];
    last_mut.post_output = px;
  }
  return self;
}

bool MetalInception::Predict(const float* input, int batch_size,
                              float* output) {
  auto& I = *impl_;

  // Fast path: no deterministic layers, run the full MPSGraph.
  if (I.det_layers.empty()) {
    int unused = 0;
    // DV_METAL_GPU_FINALIZE=1: route through "probs" tap (post dense +
    // softmax) so caller gets (B,3) probabilities directly. Default
    // path keeps "gap" which yields (B,2048) features for BnnsFinalize.
    const char* tap = I.gpu_finalize ? "probs" : "gap";
    return PredictAtTap(tap, input, batch_size, output, &unused);
  }

  // Det path: dispatch deterministic kernels in chain, then run the
  // post-graph from the last det layer's output to gap.
  @autoreleasepool {
    // 1) Allocate input buffer (user-supplied input, copied to GPU).
    const DetLayer& det0 = I.det_layers.front();
    const NSUInteger n_in =
        (NSUInteger)batch_size * det0.H_in * det0.W_in * det0.C_in;
    id<MTLBuffer> cur_buf = [I.device
        newBufferWithBytes:input
                    length:n_in * sizeof(float)
                   options:MTLResourceStorageModeShared];
    if (!cur_buf) {
      LOG(ERROR) << "MetalInception::Predict(det): input buffer alloc failed";
      return false;
    }

    // 2) Chain det stages in sequence on a single command buffer.
    id<MTLCommandBuffer> cb = [I.queue commandBuffer];
    for (size_t i = 0; i < I.det_layers.size(); ++i) {
      const DetLayer& det = I.det_layers[i];
      const NSUInteger n_dst =
          (NSUInteger)batch_size * det.H_out * det.W_out * det.C_out;
      id<MTLBuffer> dst_buf = [I.device
          newBufferWithLength:n_dst * sizeof(float)
                       options:MTLResourceStorageModeShared];
      if (!dst_buf) {
        LOG(ERROR) << "MetalInception::Predict(det): out buffer alloc failed "
                   << "for stage " << det.tap_name;
        return false;
      }
      bool ok = false;
      if (det.kind == DetLayer::kConv) {
        ConvDesc d = det.conv_desc;
        d.B = batch_size;
        d.H_in = det.H_in;
        d.W_in = det.W_in;
        d.H_out = det.H_out;
        d.W_out = det.W_out;
        if (det.use_unfolded_bn) {
          // Phase 5.5f: conv (raw, no bias, no ReLU) → bn_relu (with ReLU).
          // Use a separate intermediate buffer for the conv output, then
          // BN+ReLU writes the final dst_buf.
          id<MTLBuffer> inter_buf = [I.device
              newBufferWithLength:n_dst * sizeof(float)
                           options:MTLResourceStorageModeShared];
          if (!inter_buf) {
            LOG(ERROR) << "MetalInception::Predict(det): inter buffer alloc "
                       << "failed for stage " << det.tap_name;
            return false;
          }
          ok = I.conv_serial->Encode(cb, cur_buf, det.weights_buf,
                                      det.bias_buf, inter_buf, d);
          if (ok) {
            BnReluDesc bn_d{};
            bn_d.B = batch_size;
            bn_d.H = det.H_out;
            bn_d.W = det.W_out;
            bn_d.C = det.C_out;
            bn_d.eps = kBNEpsilon;
            bn_d.relu = true;   // post-BN ReLU
            ok = I.bn_relu->Encode(cb, inter_buf, det.bn_mean_buf,
                                     det.bn_var_buf, det.bn_beta_buf,
                                     dst_buf, bn_d);
          }
        } else {
          ok = I.conv_serial->Encode(cb, cur_buf, det.weights_buf,
                                      det.bias_buf, dst_buf, d);
        }
      } else {
        MaxPoolDesc d = det.pool_desc;
        d.B = batch_size;
        d.H_in = det.H_in;
        d.W_in = det.W_in;
        d.H_out = det.H_out;
        d.W_out = det.W_out;
        d.C = det.C_in;
        ok = I.max_pool->Encode(cb, cur_buf, dst_buf, d);
      }
      if (!ok) {
        LOG(ERROR) << "MetalInception::Predict(det): kernel encode failed "
                   << "for stage " << det.tap_name;
        return false;
      }
      cur_buf = dst_buf;  // chain
    }
    [cb commit];
    [cb waitUntilCompleted];
    if (cb.status == MTLCommandBufferStatusError) {
      LOG(ERROR) << "MetalInception::Predict(det): GPU command buffer failed: "
                 << (cb.error ? cb.error.localizedDescription.UTF8String
                              : "unknown");
      return false;
    }

    // Phase 8 / Tier 6.0 — full-network det path: bypass MPSGraph entirely
    // by chaining all 11 Inception blocks + global avg pool, then copying
    // the result to `output`. Only active when det_blocks is non-empty
    // (built when DV_METAL_SERIAL_FULL=1 + full stem chain + unfolded BN).
    if (!I.det_blocks.empty()) {
      id<MTLCommandBuffer> cb_blk = [I.queue commandBuffer];
      id<MTLBuffer> blk_in = cur_buf;
      for (size_t i = 0; i < I.det_blocks.size(); ++i) {
        if (!DispatchDetMixedBlock(cb_blk, I.conv_serial.get(),
                                    I.bn_relu.get(), I.avg_pool.get(),
                                    I.max_pool.get(), I.concat.get(),
                                    I.det_blocks[i], blk_in, batch_size)) {
          LOG(ERROR) << "MetalInception::Predict(SERIAL_FULL): block "
                     << I.det_blocks[i].tap_name << " failed";
          return false;
        }
        blk_in = I.det_blocks[i].concat_out;
      }
      // Global avg pool: (B, H, W, C) → (B, C). Last block output spatial:
      // 1×5 (after 7c on DV pileup geometry).
      const DetMixedBlock& last_blk = I.det_blocks.back();
      GlobalAvgPoolDesc gap_d{};
      gap_d.B = batch_size;
      gap_d.H_in = last_blk.H_out;
      gap_d.W_in = last_blk.W_out;
      gap_d.C = last_blk.C_out;
      if (!I.gap_pool->Encode(cb_blk, blk_in, I.gap_out_buf, gap_d)) {
        LOG(ERROR) << "MetalInception::Predict(SERIAL_FULL): gap failed";
        return false;
      }
      [cb_blk commit];
      [cb_blk waitUntilCompleted];
      if (cb_blk.status == MTLCommandBufferStatusError) {
        LOG(ERROR) << "MetalInception::Predict(det): GPU command buffer failed: "
                   << (cb_blk.error
                           ? cb_blk.error.localizedDescription.UTF8String
                           : "unknown");
        return false;
      }
      // Copy (B, feature_dim) FP32 result to `output`.
      const size_t out_bytes = (size_t)batch_size * I.feature_dim * sizeof(float);
      std::memcpy(output, [I.gap_out_buf contents], out_bytes);
      return true;
    }

    // 3) Run post-graph (last det layer output → gap) using cur_buf as
    // the placeholder.
    const DetLayer& last = I.det_layers.back();
    NSNumber* bs_key = @(batch_size);
    MPSGraphExecutable* post_exe = I.post_exec_cache[bs_key];
    if (!post_exe) {
      MPSShape* in_shape = @[@(batch_size),
                              @(last.H_out),
                              @(last.W_out),
                              @(last.C_out)];
      MPSGraphShapedType* in_st =
          [[MPSGraphShapedType alloc] initWithShape:in_shape
                                            dataType:MPSDataTypeFloat32];
      NSDictionary<MPSGraphTensor*, MPSGraphShapedType*>* feeds_shape =
          @{last.post_input: in_st};
      post_exe = [last.post_graph
          compileWithDevice:[MPSGraphDevice deviceWithMTLDevice:I.device]
                      feeds:feeds_shape
              targetTensors:@[last.post_output]
           targetOperations:nil
      compilationDescriptor:I.compileDesc];
      if (!post_exe) {
        LOG(ERROR) << "MetalInception::Predict(det): post-graph compile failed";
        return false;
      }
      I.post_exec_cache[bs_key] = post_exe;
    }

    MPSGraphTensorData* in_td = [[MPSGraphTensorData alloc]
        initWithMTLBuffer:cur_buf
                    shape:@[@(batch_size),
                             @(last.H_out),
                             @(last.W_out),
                             @(last.C_out)]
                 dataType:MPSDataTypeFloat32];

    MPSGraphExecutableExecutionDescriptor* runDesc =
        [MPSGraphExecutableExecutionDescriptor new];
    runDesc.waitUntilCompleted = YES;
    NSArray<MPSGraphTensorData*>* outs =
        [post_exe runWithMTLCommandQueue:I.queue
                            inputsArray:@[in_td]
                           resultsArray:nil
                    executionDescriptor:runDesc];
    if (!outs || outs.count != 1) {
      LOG(ERROR) << "MetalInception::Predict(det): post-graph run produced "
                 << (outs ? outs.count : 0) << " results";
      return false;
    }
    [outs[0].mpsndarray readBytes:output strideBytes:nil];
  }
  return true;
}

bool MetalInception::PredictAtTap(const std::string& tap_name,
                                   const float* input, int batch_size,
                                   float* output,
                                   int* out_total_elems_per_image) {
  if (!input || !output || batch_size <= 0) {
    LOG(ERROR) << "MetalInception::PredictAtTap: bad args";
    return false;
  }
  auto& I = *impl_;

  @autoreleasepool {
    NSString* tap_ns = [NSString stringWithUTF8String:tap_name.c_str()];
    MPSGraphTensor* tap = I.taps[tap_ns];
    if (!tap) {
      LOG(ERROR) << "MetalInception: unknown tap '" << tap_name << "'";
      return false;
    }

    // Compile (and cache) an executable for this specific tap with
    // optimizationLevel=Level0 — only path that gives correct FP32
    // output (Phase 5.5a investigation).
    MPSGraphExecutable* exe = I.execCache[tap_ns];
    if (!exe) {
      MPSShape* in_shape =
          @[@(batch_size), @(I.input_height), @(I.input_width), @(I.input_channels)];
      MPSGraphShapedType* in_st =
          [[MPSGraphShapedType alloc] initWithShape:in_shape
                                            dataType:MPSDataTypeFloat32];
      NSDictionary<MPSGraphTensor*, MPSGraphShapedType*>* feeds_shape =
          @{I.input: in_st};
      exe = [I.graph compileWithDevice:[MPSGraphDevice deviceWithMTLDevice:I.device]
                                  feeds:feeds_shape
                          targetTensors:@[tap]
                       targetOperations:nil
                  compilationDescriptor:I.compileDesc];
      if (!exe) {
        LOG(ERROR) << "MetalInception::PredictAtTap: compile failed for "
                   << tap_name;
        return false;
      }
      I.execCache[tap_ns] = exe;
    }

    // Wrap input as MPSGraphTensorData.
    const NSUInteger n_in = (NSUInteger)batch_size *
        I.input_height * I.input_width * I.input_channels;
    NSData* in_data = [NSData dataWithBytes:input
                                     length:n_in * sizeof(float)];
    MPSGraphTensorData* in_td = [[MPSGraphTensorData alloc]
        initWithDevice:[MPSGraphDevice deviceWithMTLDevice:I.device]
                  data:in_data
                 shape:@[@(batch_size), @(I.input_height), @(I.input_width),
                          @(I.input_channels)]
              dataType:MPSDataTypeFloat32];

    MPSGraphExecutableExecutionDescriptor* runDesc =
        [MPSGraphExecutableExecutionDescriptor new];
    runDesc.waitUntilCompleted = YES;
    NSArray<MPSGraphTensorData*>* outs =
        [exe runWithMTLCommandQueue:I.queue
                        inputsArray:@[in_td]
                       resultsArray:nil
                executionDescriptor:runDesc];
    if (!outs || outs.count != 1) {
      LOG(ERROR) << "MetalInception::PredictAtTap: run produced "
                 << (outs ? outs.count : 0) << " results (expected 1)";
      return false;
    }
    MPSGraphTensorData* out_td = outs[0];
    NSArray<NSNumber*>* shape = out_td.shape;
    NSUInteger total = 1;
    for (NSNumber* d in shape) total *= [d unsignedIntegerValue];
    if (out_total_elems_per_image && batch_size > 0) {
      *out_total_elems_per_image =
          static_cast<int>(total / (NSUInteger)batch_size);
    }
    [out_td.mpsndarray readBytes:output strideBytes:nil];
  }
  return true;
}

}  // namespace deepvariant
