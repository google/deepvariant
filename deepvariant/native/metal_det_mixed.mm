// Phase 8 / Tier 6.0 — Deterministic Inception block dispatch for all
// 11 Mixed_X blocks (5b, 5c, 5d, 6a, 6b-6e, 7a, 7b-7c).
//
// Each block builder loads .dvw weights for its (conv_n, bn_n) tuples
// and allocates intermediate + output MTLBuffers. The unified
// DispatchDetMixedBlock encoder dispatches all branch ops + concat
// onto a single MTLCommandBuffer, supporting:
//   - Sequential branches (Mixed_5b/5c/5d, 6b-6e, etc.)
//   - Pool-only branches (Mixed_6a/7a max-pool branch)
//   - Avg-pool prepended branches (Mixed_5x/6b-e/7b-c pool branch)
//   - Split branches (Mixed_7b/7c b3a, b3b — trunk + 1×3 + 3×1 split)

#include "deepvariant/native/metal_det_mixed.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <vector>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr float kBNEpsilon = 1e-3f;

// When true, BuildBranchOp folds BN into conv weights at build time
// (matches the baseline MPSGraph CBR path that produces 160 FM vs Docker
// at chr20-full scale). When false, builds raw conv + separate BN+ReLU
// (the UNFOLDED path that produces 8837 FM regression at chr20-full scale
// — provided only as a research toggle via DV_METAL_DET_UNFOLDED=1).
bool g_det_use_folded = true;

std::string DetAttr(int n, const char* attr) {
  return "layer_with_weights-" + std::to_string(n) +
         "/" + attr + "/.ATTRIBUTES/VARIABLE_VALUE";
}

// Local copy of the FoldConvBn pattern from metal_inference.mm. Folds
// (Conv HWIO + BN gamma=1, beta, mean, var, eps) into a single
// (W' HWIO, b') pair where:
//   scale[o] = 1 / sqrt(var[o] + eps)
//   W'[h,w,i,o] = W[h,w,i,o] * scale[o]
//   b'[o] = beta[o] - mean[o] * scale[o]
struct DetFusedConv {
  std::vector<float> weights_hwio;
  std::vector<float> bias;
  int O = 0, I = 0, H = 0, W = 0;
};
DetFusedConv DetFoldConvBn(const DvwWeights& dvw, int conv_n, int bn_n) {
  const auto* k = dvw.Get(DetAttr(conv_n, "kernel"));
  const auto* beta = dvw.Get(DetAttr(bn_n, "beta"));
  const auto* mean = dvw.Get(DetAttr(bn_n, "moving_mean"));
  const auto* var = dvw.Get(DetAttr(bn_n, "moving_variance"));
  if (!k || !beta || !mean || !var || k->shape.size() != 4u) return {};
  const int Hk = k->shape[0], Wk = k->shape[1];
  const int Ik = k->shape[2], Ok = k->shape[3];
  DetFusedConv out;
  out.H = Hk; out.W = Wk; out.I = Ik; out.O = Ok;
  std::vector<float> scale(Ok), offset(Ok);
  for (int o = 0; o < Ok; ++o) {
    scale[o] = 1.0f / std::sqrt(var->data[o] + kBNEpsilon);
    offset[o] = beta->data[o] - mean->data[o] * scale[o];
  }
  out.bias = std::move(offset);
  out.weights_hwio.resize((size_t)Hk * Wk * Ik * Ok);
  for (size_t h = 0; h < (size_t)Hk; ++h) {
    for (size_t w = 0; w < (size_t)Wk; ++w) {
      for (size_t i = 0; i < (size_t)Ik; ++i) {
        for (size_t o = 0; o < (size_t)Ok; ++o) {
          const size_t idx = ((h * Wk + w) * Ik + i) * Ok + o;
          out.weights_hwio[idx] = k->data[idx] * scale[o];
        }
      }
    }
  }
  return out;
}

id<MTLBuffer> NewBuf(id<MTLDevice> dev, size_t bytes) {
  return [dev newBufferWithLength:bytes
                          options:MTLResourceStorageModeShared];
}

// Build one CBR (conv + BN + ReLU) op. By default uses FOLDED weights
// (W' = W*scale, bias = offset, fused ReLU in conv) — bit-equivalent
// to the baseline MPSGraph CBR path that achieves 100 % FILTER parity
// on chr20:10M-10.1M and 160 FM on chr20 full HG003.
//
// The unfolded path (raw conv + separate MetalBnRelu) is bit-different
// at scale (8837 FM regression confirmed in chr20 full HG003 testing —
// same magnitude as Probe C2 unfolded MPSGraph) and is provided only
// as a research toggle via g_det_use_folded = false.
bool BuildBranchOp(id<MTLDevice> device, const DvwWeights& dvw,
                   int conv_n, int bn_n,
                   int H_in, int W_in, int C_in,
                   int H_out, int W_out, int stride_h, int stride_w,
                   bool same_padding, int max_B,
                   DetBranchOp* op) {
  if (g_det_use_folded) {
    // ── FOLDED path ── conv emits final activation with bias + ReLU
    // applied in the kernel. No mean/var/beta needed at runtime.
    DetFusedConv fc = DetFoldConvBn(dvw, conv_n, bn_n);
    if (fc.weights_hwio.empty()) {
      LOG(ERROR) << "BuildBranchOp: FoldConvBn failed for conv=" << conv_n
                 << " bn=" << bn_n;
      return false;
    }
    if (fc.I != C_in) {
      LOG(ERROR) << "BuildBranchOp(folded): weight C_in=" << fc.I
                 << " mismatch geom C_in=" << C_in
                 << " (conv=" << conv_n << ")";
      return false;
    }
    op->conv.B = max_B;
    op->conv.H_in = H_in;     op->conv.W_in = W_in;     op->conv.C_in = C_in;
    op->conv.H_out = H_out;   op->conv.W_out = W_out;   op->conv.C_out = fc.O;
    op->conv.Kh = fc.H;       op->conv.Kw = fc.W;
    op->conv.stride_h = stride_h; op->conv.stride_w = stride_w;
    op->conv.pad_h = same_padding ? (fc.H - 1) / 2 : 0;
    op->conv.pad_w = same_padding ? (fc.W - 1) / 2 : 0;
    op->conv.relu = true;     // fused ReLU
    op->w =
        [device newBufferWithBytes:fc.weights_hwio.data()
                            length:fc.weights_hwio.size() * sizeof(float)
                           options:MTLResourceStorageModeShared];
    op->bias =
        [device newBufferWithBytes:fc.bias.data()
                            length:fc.bias.size() * sizeof(float)
                           options:MTLResourceStorageModeShared];
    op->mean = nil;     // unused in folded path
    op->var = nil;
    op->beta = nil;
    op->raw_buf = nil;  // no separate BN intermediate
    const size_t act_bytes = (size_t)max_B * H_out * W_out * fc.O * sizeof(float);
    op->out_buf = NewBuf(device, act_bytes);
    if (!op->w || !op->bias || !op->out_buf) {
      LOG(ERROR) << "BuildBranchOp(folded): alloc failed for conv=" << conv_n;
      return false;
    }
    op->out_H = H_out; op->out_W = W_out; op->out_C = fc.O;
    return true;
  }

  // ── UNFOLDED path (research toggle) ─── raw conv + separate BN+ReLU.
  const auto* k = dvw.Get(DetAttr(conv_n, "kernel"));
  const auto* beta = dvw.Get(DetAttr(bn_n, "beta"));
  const auto* mean = dvw.Get(DetAttr(bn_n, "moving_mean"));
  const auto* var = dvw.Get(DetAttr(bn_n, "moving_variance"));
  if (!k || !beta || !mean || !var || k->shape.size() != 4u) {
    LOG(ERROR) << "BuildBranchOp(unfolded): missing weights for conv=" << conv_n
               << " bn=" << bn_n;
    return false;
  }
  const int Hk = k->shape[0], Wk = k->shape[1];
  const int Ik = k->shape[2], Ok = k->shape[3];
  if (Ik != C_in) return false;
  op->conv.B = max_B;
  op->conv.H_in = H_in; op->conv.W_in = W_in; op->conv.C_in = C_in;
  op->conv.H_out = H_out; op->conv.W_out = W_out; op->conv.C_out = Ok;
  op->conv.Kh = Hk; op->conv.Kw = Wk;
  op->conv.stride_h = stride_h; op->conv.stride_w = stride_w;
  op->conv.pad_h = same_padding ? (Hk - 1) / 2 : 0;
  op->conv.pad_w = same_padding ? (Wk - 1) / 2 : 0;
  op->conv.relu = false;
  op->w = [device newBufferWithBytes:k->data length:k->n_bytes
                             options:MTLResourceStorageModeShared];
  std::vector<float> zero_bias(Ok, 0.0f);
  op->bias = [device newBufferWithBytes:zero_bias.data()
                                 length:Ok * sizeof(float)
                                options:MTLResourceStorageModeShared];
  op->mean = [device newBufferWithBytes:mean->data length:Ok * sizeof(float)
                                options:MTLResourceStorageModeShared];
  op->var = [device newBufferWithBytes:var->data length:Ok * sizeof(float)
                               options:MTLResourceStorageModeShared];
  op->beta = [device newBufferWithBytes:beta->data length:Ok * sizeof(float)
                                options:MTLResourceStorageModeShared];
  const size_t act_bytes = (size_t)max_B * H_out * W_out * Ok * sizeof(float);
  op->raw_buf = NewBuf(device, act_bytes);
  op->out_buf = NewBuf(device, act_bytes);
  if (!op->w || !op->bias || !op->mean || !op->var || !op->beta ||
      !op->raw_buf || !op->out_buf) return false;
  op->out_H = H_out; op->out_W = W_out; op->out_C = Ok;
  return true;
}

// Build an InceptionA-style block (Mixed_5b/5c/5d): 4 branches, all
// SAME padding, no spatial change.
//   br0: 1×1                         -> C_b1
//   br1: 1×1 -> 5×5                  -> C_b5
//   br2: 1×1 -> 3×3 -> 3×3           -> C_b3
//   br3: avg_pool 3×3 -> 1×1         -> C_bp
bool BuildInceptionA(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     int b1_conv, int b1_bn,
                     int b5a_conv, int b5a_bn, int b5b_conv, int b5b_bn,
                     int b3a_conv, int b3a_bn,
                     int b3b_conv, int b3b_bn,
                     int b3c_conv, int b3c_bn,
                     int bp_conv, int bp_bn,
                     const std::string& tap, DetMixedBlock* block) {
  block->tap_name = tap;
  block->B = max_B;
  block->H_in = H_in;
  block->W_in = W_in;
  block->C_in = C_in;
  block->H_out = H_in;
  block->W_out = W_in;
  block->branches.clear();
  block->branches.resize(4);

  // br0: 1×1
  block->branches[0].ops.resize(1);
  if (!BuildBranchOp(device, dvw, b1_conv, b1_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[0].ops[0])) return false;

  // br1: 1×1 -> 5×5
  block->branches[1].ops.resize(2);
  if (!BuildBranchOp(device, dvw, b5a_conv, b5a_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[1].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, b5b_conv, b5b_bn,
                     H_in, W_in, block->branches[1].ops[0].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &block->branches[1].ops[1])) return false;

  // br2: 1×1 -> 3×3 -> 3×3
  block->branches[2].ops.resize(3);
  if (!BuildBranchOp(device, dvw, b3a_conv, b3a_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[2].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, b3b_conv, b3b_bn,
                     H_in, W_in, block->branches[2].ops[0].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &block->branches[2].ops[1])) return false;
  if (!BuildBranchOp(device, dvw, b3c_conv, b3c_bn,
                     H_in, W_in, block->branches[2].ops[1].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &block->branches[2].ops[2])) return false;

  // br3: avg_pool 3×3 -> 1×1
  {
    DetBranch& br = block->branches[3];
    br.has_avg_pool_pre = true;
    br.avg_pool.B = max_B;
    br.avg_pool.H_in = H_in; br.avg_pool.W_in = W_in;
    br.avg_pool.C = C_in;
    br.avg_pool.H_out = H_in; br.avg_pool.W_out = W_in;
    br.avg_pool.Kh = 3; br.avg_pool.Kw = 3;
    br.avg_pool.stride_h = 1; br.avg_pool.stride_w = 1;
    br.avg_pool.pad_h = 1; br.avg_pool.pad_w = 1;
    br.avg_pool.exclude_pad = true;
    const size_t pool_bytes = (size_t)max_B * H_in * W_in * C_in * sizeof(float);
    br.pool_out = NewBuf(device, pool_bytes);
    if (!br.pool_out) return false;
    br.ops.resize(1);
    if (!BuildBranchOp(device, dvw, bp_conv, bp_bn,
                       H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                       &br.ops[0])) return false;
  }

  block->C_out =
      block->branches[0].ops.back().out_C +
      block->branches[1].ops.back().out_C +
      block->branches[2].ops.back().out_C +
      block->branches[3].ops.back().out_C;
  const size_t concat_bytes =
      (size_t)max_B * block->H_out * block->W_out * block->C_out * sizeof(float);
  block->concat_out = NewBuf(device, concat_bytes);
  if (!block->concat_out) return false;
  return true;
}

// Build an InceptionB-style block (Mixed_6b/6c/6d/6e): 4 branches,
// SAME padding, no spatial change. Asymmetric 7×7 factorisation.
//   br0: 1×1                              -> C_b1
//   br1: 1×1 -> 1×7 -> 7×1                -> C_b7a
//   br2: 1×1 -> 7×1 -> 1×7 -> 7×1 -> 1×7  -> C_b7b
//   br3: avg_pool 3×3 -> 1×1              -> C_bp
bool BuildInceptionB(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     int b1_conv, int b1_bn,
                     int b7a0_conv, int b7a0_bn,
                     int b7a1_conv, int b7a1_bn,
                     int b7a2_conv, int b7a2_bn,
                     int b7b0_conv, int b7b0_bn,
                     int b7b1_conv, int b7b1_bn,
                     int b7b2_conv, int b7b2_bn,
                     int b7b3_conv, int b7b3_bn,
                     int b7b4_conv, int b7b4_bn,
                     int bp_conv, int bp_bn,
                     const std::string& tap, DetMixedBlock* block) {
  block->tap_name = tap;
  block->B = max_B;
  block->H_in = H_in; block->W_in = W_in; block->C_in = C_in;
  block->H_out = H_in; block->W_out = W_in;
  block->branches.clear();
  block->branches.resize(4);

  // br0: 1×1
  block->branches[0].ops.resize(1);
  if (!BuildBranchOp(device, dvw, b1_conv, b1_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[0].ops[0])) return false;

  // br1: 1×1 -> 1×7 -> 7×1
  block->branches[1].ops.resize(3);
  if (!BuildBranchOp(device, dvw, b7a0_conv, b7a0_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[1].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, b7a1_conv, b7a1_bn,
                     H_in, W_in, block->branches[1].ops[0].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &block->branches[1].ops[1])) return false;
  if (!BuildBranchOp(device, dvw, b7a2_conv, b7a2_bn,
                     H_in, W_in, block->branches[1].ops[1].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &block->branches[1].ops[2])) return false;

  // br2: 1×1 -> 7×1 -> 1×7 -> 7×1 -> 1×7
  block->branches[2].ops.resize(5);
  if (!BuildBranchOp(device, dvw, b7b0_conv, b7b0_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &block->branches[2].ops[0])) return false;
  for (int i = 1; i < 5; ++i) {
    int conv_n = (i == 1) ? b7b1_conv :
                 (i == 2) ? b7b2_conv :
                 (i == 3) ? b7b3_conv : b7b4_conv;
    int bn_n = (i == 1) ? b7b1_bn :
               (i == 2) ? b7b2_bn :
               (i == 3) ? b7b3_bn : b7b4_bn;
    if (!BuildBranchOp(device, dvw, conv_n, bn_n,
                       H_in, W_in, block->branches[2].ops[i-1].out_C,
                       H_in, W_in, 1, 1, true, max_B,
                       &block->branches[2].ops[i])) return false;
  }

  // br3: avg_pool 3×3 -> 1×1
  {
    DetBranch& br = block->branches[3];
    br.has_avg_pool_pre = true;
    br.avg_pool.B = max_B;
    br.avg_pool.H_in = H_in; br.avg_pool.W_in = W_in;
    br.avg_pool.C = C_in;
    br.avg_pool.H_out = H_in; br.avg_pool.W_out = W_in;
    br.avg_pool.Kh = 3; br.avg_pool.Kw = 3;
    br.avg_pool.stride_h = 1; br.avg_pool.stride_w = 1;
    br.avg_pool.pad_h = 1; br.avg_pool.pad_w = 1;
    br.avg_pool.exclude_pad = true;
    const size_t pool_bytes = (size_t)max_B * H_in * W_in * C_in * sizeof(float);
    br.pool_out = NewBuf(device, pool_bytes);
    if (!br.pool_out) return false;
    br.ops.resize(1);
    if (!BuildBranchOp(device, dvw, bp_conv, bp_bn,
                       H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                       &br.ops[0])) return false;
  }

  block->C_out =
      block->branches[0].ops.back().out_C +
      block->branches[1].ops.back().out_C +
      block->branches[2].ops.back().out_C +
      block->branches[3].ops.back().out_C;
  const size_t concat_bytes =
      (size_t)max_B * block->H_out * block->W_out * block->C_out * sizeof(float);
  block->concat_out = NewBuf(device, concat_bytes);
  return block->concat_out != nil;
}

}  // namespace

// =============================================================
//  Block builders
// =============================================================

bool BuildDetMixed5b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  // Indices from metal_inference.mm Mixed_5b().
  return BuildInceptionA(device, dvw, max_B, H_in, W_in, C_in,
      /*b1*/ 16, 20,
      /*b5a*/ 12, 14, /*b5b*/ 17, 21,
      /*b3a*/ 10, 11, /*b3b*/ 13, 15, /*b3c*/ 18, 22,
      /*bp*/ 19, 23, "5b", out);
}

bool BuildDetMixed5c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionA(device, dvw, max_B, H_in, W_in, C_in,
      30, 34,
      26, 28, 31, 35,
      24, 25, 27, 29, 32, 36,
      33, 37, "5c", out);
}

bool BuildDetMixed5d(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionA(device, dvw, max_B, H_in, W_in, C_in,
      44, 48,
      40, 42, 45, 49,
      38, 39, 41, 43, 46, 50,
      47, 51, "5d", out);
}

// Mixed_6a (Reduction-A): 3 branches, stride-2 VALID.
//   br0: 3×3 stride-2 VALID 288→384       -> C=384
//   br1: 1×1 SAME -> 3×3 SAME -> 3×3 stride-2 VALID
//   br2: max_pool 3×3 stride-2 VALID      -> C_in (passes through)
// Output spatial: H_in→ceil((H_in - 3 + 1) / 2) (VALID).
bool BuildDetMixed6a(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  out->tap_name = "6a";
  out->B = max_B;
  out->H_in = H_in; out->W_in = W_in; out->C_in = C_in;
  // VALID 3×3 stride 2: H_out = floor((H_in - 3)/2) + 1
  out->H_out = (H_in - 3) / 2 + 1;
  out->W_out = (W_in - 3) / 2 + 1;
  out->branches.clear();
  out->branches.resize(3);

  // br0: 3×3 stride-2 VALID
  out->branches[0].ops.resize(1);
  if (!BuildBranchOp(device, dvw, 56, 58,
                     H_in, W_in, C_in, out->H_out, out->W_out, 2, 2, false, max_B,
                     &out->branches[0].ops[0])) return false;

  // br1: 1×1 SAME -> 3×3 SAME -> 3×3 stride-2 VALID
  out->branches[1].ops.resize(3);
  if (!BuildBranchOp(device, dvw, 52, 53,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &out->branches[1].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, 54, 55,
                     H_in, W_in, out->branches[1].ops[0].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &out->branches[1].ops[1])) return false;
  if (!BuildBranchOp(device, dvw, 57, 59,
                     H_in, W_in, out->branches[1].ops[1].out_C,
                     out->H_out, out->W_out, 2, 2, false, max_B,
                     &out->branches[1].ops[2])) return false;

  // br2: max-pool 3×3 stride-2 VALID — pool-only, no convs.
  {
    DetBranch& br = out->branches[2];
    br.pool_only = true;
    br.has_max_pool_pre = false;  // pool IS the branch op; not a "pre"
    br.max_pool.B = max_B;
    br.max_pool.H_in = H_in; br.max_pool.W_in = W_in; br.max_pool.C = C_in;
    br.max_pool.H_out = out->H_out; br.max_pool.W_out = out->W_out;
    br.max_pool.Kh = 3; br.max_pool.Kw = 3;
    br.max_pool.stride_h = 2; br.max_pool.stride_w = 2;
    br.max_pool.pad_h = 0; br.max_pool.pad_w = 0;
    const size_t pool_bytes =
        (size_t)max_B * out->H_out * out->W_out * C_in * sizeof(float);
    br.pool_out = NewBuf(device, pool_bytes);
    if (!br.pool_out) return false;
    br.split_out_C = C_in;  // re-used field as "branch output channel count"
  }

  out->C_out = out->branches[0].ops.back().out_C +
               out->branches[1].ops.back().out_C +
               C_in;
  const size_t concat_bytes =
      (size_t)max_B * out->H_out * out->W_out * out->C_out * sizeof(float);
  out->concat_out = NewBuf(device, concat_bytes);
  return out->concat_out != nil;
}

bool BuildDetMixed6b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionB(device, dvw, max_B, H_in, W_in, C_in,
      /*b1*/ 72, 76,
      /*b7a*/ 64, 66, 68, 70, 73, 77,
      /*b7b*/ 60, 61, 62, 63, 65, 67, 69, 71, 74, 78,
      /*bp*/ 75, 79, "6b", out);
}

bool BuildDetMixed6c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionB(device, dvw, max_B, H_in, W_in, C_in,
      92, 96,
      84, 86, 88, 90, 93, 97,
      80, 81, 82, 83, 85, 87, 89, 91, 94, 98,
      95, 99, "6c", out);
}

bool BuildDetMixed6d(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionB(device, dvw, max_B, H_in, W_in, C_in,
      112, 116,
      104, 106, 108, 110, 113, 117,
      100, 101, 102, 103, 105, 107, 109, 111, 114, 118,
      115, 119, "6d", out);
}

bool BuildDetMixed6e(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionB(device, dvw, max_B, H_in, W_in, C_in,
      132, 136,
      124, 126, 128, 130, 133, 137,
      120, 121, 122, 123, 125, 127, 129, 131, 134, 138,
      135, 139, "6e", out);
}

// Mixed_7a (Reduction-B): 3 branches, stride-2 VALID at the end.
//   br0: 1×1 SAME -> 3×3 stride-2 VALID
//   br1: 1×1 SAME -> 1×7 SAME -> 7×1 SAME -> 3×3 stride-2 VALID
//   br2: max_pool 3×3 stride-2 VALID
bool BuildDetMixed7a(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  out->tap_name = "7a";
  out->B = max_B;
  out->H_in = H_in; out->W_in = W_in; out->C_in = C_in;
  out->H_out = (H_in - 3) / 2 + 1;
  out->W_out = (W_in - 3) / 2 + 1;
  out->branches.clear();
  out->branches.resize(3);

  // br0: 1×1 SAME -> 3×3 stride-2 VALID
  out->branches[0].ops.resize(2);
  if (!BuildBranchOp(device, dvw, 144, 146,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &out->branches[0].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, 148, 150,
                     H_in, W_in, out->branches[0].ops[0].out_C,
                     out->H_out, out->W_out, 2, 2, false, max_B,
                     &out->branches[0].ops[1])) return false;

  // br1: 1×1 -> 1×7 -> 7×1 -> 3×3 stride-2 VALID
  out->branches[1].ops.resize(4);
  if (!BuildBranchOp(device, dvw, 140, 141,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &out->branches[1].ops[0])) return false;
  if (!BuildBranchOp(device, dvw, 142, 143,
                     H_in, W_in, out->branches[1].ops[0].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &out->branches[1].ops[1])) return false;
  if (!BuildBranchOp(device, dvw, 145, 147,
                     H_in, W_in, out->branches[1].ops[1].out_C,
                     H_in, W_in, 1, 1, true, max_B,
                     &out->branches[1].ops[2])) return false;
  if (!BuildBranchOp(device, dvw, 149, 151,
                     H_in, W_in, out->branches[1].ops[2].out_C,
                     out->H_out, out->W_out, 2, 2, false, max_B,
                     &out->branches[1].ops[3])) return false;

  // br2: max-pool 3×3 stride-2 VALID
  {
    DetBranch& br = out->branches[2];
    br.pool_only = true;
    br.max_pool.B = max_B;
    br.max_pool.H_in = H_in; br.max_pool.W_in = W_in; br.max_pool.C = C_in;
    br.max_pool.H_out = out->H_out; br.max_pool.W_out = out->W_out;
    br.max_pool.Kh = 3; br.max_pool.Kw = 3;
    br.max_pool.stride_h = 2; br.max_pool.stride_w = 2;
    br.max_pool.pad_h = 0; br.max_pool.pad_w = 0;
    const size_t pool_bytes =
        (size_t)max_B * out->H_out * out->W_out * C_in * sizeof(float);
    br.pool_out = NewBuf(device, pool_bytes);
    if (!br.pool_out) return false;
    br.split_out_C = C_in;
  }

  out->C_out = out->branches[0].ops.back().out_C +
               out->branches[1].ops.back().out_C +
               C_in;
  const size_t concat_bytes =
      (size_t)max_B * out->H_out * out->W_out * out->C_out * sizeof(float);
  out->concat_out = NewBuf(device, concat_bytes);
  return out->concat_out != nil;
}

// Mixed_7b/7c (InceptionC): 4 branches, two with split outputs.
//   br0: 1×1                                                       -> 320 ch
//   br1 (split): 1×1 -> {1×3, 3×1}                                  -> 384+384 = 768 ch
//   br2 (split): 1×1 -> 3×3 -> {1×3, 3×1}                          -> 384+384 = 768 ch
//   br3: avg_pool 3×3 -> 1×1                                       -> 192 ch
// Note ops layout for split branches: trunk_size = N (sequential prefix),
// then ops[N] and ops[N+1] are PARALLEL on trunk's last out_buf.
static bool BuildInceptionC(id<MTLDevice> device, const DvwWeights& dvw,
                            int max_B, int H_in, int W_in, int C_in,
                            int b1_conv, int b1_bn,
                            int b3a_trunk_conv, int b3a_trunk_bn,
                            int b3a_1x3_conv, int b3a_1x3_bn,
                            int b3a_3x1_conv, int b3a_3x1_bn,
                            int b3b_t0_conv, int b3b_t0_bn,
                            int b3b_t1_conv, int b3b_t1_bn,
                            int b3b_1x3_conv, int b3b_1x3_bn,
                            int b3b_3x1_conv, int b3b_3x1_bn,
                            int bp_conv, int bp_bn,
                            const std::string& tap, DetMixedBlock* out) {
  out->tap_name = tap;
  out->B = max_B;
  out->H_in = H_in; out->W_in = W_in; out->C_in = C_in;
  out->H_out = H_in; out->W_out = W_in;
  out->branches.clear();
  out->branches.resize(4);

  // br0: 1×1
  out->branches[0].ops.resize(1);
  if (!BuildBranchOp(device, dvw, b1_conv, b1_bn,
                     H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                     &out->branches[0].ops[0])) return false;

  // br1: split 1×1 -> {1×3, 3×1}
  {
    DetBranch& br = out->branches[1];
    br.is_split = true;
    br.trunk_size = 1;
    br.ops.resize(3);
    if (!BuildBranchOp(device, dvw, b3a_trunk_conv, b3a_trunk_bn,
                       H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                       &br.ops[0])) return false;
    if (!BuildBranchOp(device, dvw, b3a_1x3_conv, b3a_1x3_bn,
                       H_in, W_in, br.ops[0].out_C,
                       H_in, W_in, 1, 1, true, max_B, &br.ops[1])) return false;
    if (!BuildBranchOp(device, dvw, b3a_3x1_conv, b3a_3x1_bn,
                       H_in, W_in, br.ops[0].out_C,
                       H_in, W_in, 1, 1, true, max_B, &br.ops[2])) return false;
    br.split_out_C = br.ops[1].out_C + br.ops[2].out_C;
    const size_t scbytes =
        (size_t)max_B * H_in * W_in * br.split_out_C * sizeof(float);
    br.split_concat_out = NewBuf(device, scbytes);
    if (!br.split_concat_out) return false;
  }

  // br2: split 1×1 -> 3×3 -> {1×3, 3×1}
  {
    DetBranch& br = out->branches[2];
    br.is_split = true;
    br.trunk_size = 2;
    br.ops.resize(4);
    if (!BuildBranchOp(device, dvw, b3b_t0_conv, b3b_t0_bn,
                       H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                       &br.ops[0])) return false;
    if (!BuildBranchOp(device, dvw, b3b_t1_conv, b3b_t1_bn,
                       H_in, W_in, br.ops[0].out_C,
                       H_in, W_in, 1, 1, true, max_B, &br.ops[1])) return false;
    if (!BuildBranchOp(device, dvw, b3b_1x3_conv, b3b_1x3_bn,
                       H_in, W_in, br.ops[1].out_C,
                       H_in, W_in, 1, 1, true, max_B, &br.ops[2])) return false;
    if (!BuildBranchOp(device, dvw, b3b_3x1_conv, b3b_3x1_bn,
                       H_in, W_in, br.ops[1].out_C,
                       H_in, W_in, 1, 1, true, max_B, &br.ops[3])) return false;
    br.split_out_C = br.ops[2].out_C + br.ops[3].out_C;
    const size_t scbytes =
        (size_t)max_B * H_in * W_in * br.split_out_C * sizeof(float);
    br.split_concat_out = NewBuf(device, scbytes);
    if (!br.split_concat_out) return false;
  }

  // br3: avg_pool 3×3 -> 1×1
  {
    DetBranch& br = out->branches[3];
    br.has_avg_pool_pre = true;
    br.avg_pool.B = max_B;
    br.avg_pool.H_in = H_in; br.avg_pool.W_in = W_in;
    br.avg_pool.C = C_in;
    br.avg_pool.H_out = H_in; br.avg_pool.W_out = W_in;
    br.avg_pool.Kh = 3; br.avg_pool.Kw = 3;
    br.avg_pool.stride_h = 1; br.avg_pool.stride_w = 1;
    br.avg_pool.pad_h = 1; br.avg_pool.pad_w = 1;
    br.avg_pool.exclude_pad = true;
    const size_t pool_bytes = (size_t)max_B * H_in * W_in * C_in * sizeof(float);
    br.pool_out = NewBuf(device, pool_bytes);
    if (!br.pool_out) return false;
    br.ops.resize(1);
    if (!BuildBranchOp(device, dvw, bp_conv, bp_bn,
                       H_in, W_in, C_in, H_in, W_in, 1, 1, true, max_B,
                       &br.ops[0])) return false;
  }

  out->C_out =
      out->branches[0].ops.back().out_C +
      out->branches[1].split_out_C +
      out->branches[2].split_out_C +
      out->branches[3].ops.back().out_C;
  const size_t concat_bytes =
      (size_t)max_B * out->H_out * out->W_out * out->C_out * sizeof(float);
  out->concat_out = NewBuf(device, concat_bytes);
  return out->concat_out != nil;
}

bool BuildDetMixed7b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionC(device, dvw, max_B, H_in, W_in, C_in,
      /*b1*/        162, 168,
      /*b3a trunk*/ 154, 156,
      /*b3a 1×3*/   158, 163,
      /*b3a 3×1*/   159, 164,
      /*b3b t0*/    152, 153,
      /*b3b t1*/    155, 157,
      /*b3b 1×3*/   160, 165,
      /*b3b 3×1*/   161, 166,
      /*bp*/        167, 169, "7b", out);
}

bool BuildDetMixed7c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out) {
  return BuildInceptionC(device, dvw, max_B, H_in, W_in, C_in,
      180, 186,
      172, 174, 176, 181, 177, 182,
      170, 171, 173, 175, 178, 183, 179, 184,
      185, 187, "7c", out);
}

// =============================================================
//  Dispatch
// =============================================================

bool DispatchDetMixedBlock(id<MTLCommandBuffer> cb,
                           MetalConvSerial* conv_serial,
                           MetalBnRelu* bn_relu,
                           MetalAvgPool* avg_pool,
                           MetalMaxPool* max_pool,
                           MetalConcat* concat,
                           const DetMixedBlock& block,
                           id<MTLBuffer> input_buf,
                           int batch_size) {
  // bn_relu may be nil in folded mode (BN baked into conv weights).
  // max_pool may be nil for blocks without max-pool branch.
  if (!cb || !conv_serial || !avg_pool || !concat || !input_buf) {
    LOG(ERROR) << "DispatchDetMixedBlock: nil arg";
    return false;
  }
  if (batch_size <= 0 || batch_size > block.B) {
    LOG(ERROR) << "DispatchDetMixedBlock: bad batch_size=" << batch_size;
    return false;
  }

  // Per-branch dispatch; collect branch output buffer pointers + channel counts.
  id<MTLBuffer> br_outs[4] = {nil, nil, nil, nil};
  int br_c[4] = {0, 0, 0, 0};

  for (size_t bi = 0; bi < block.branches.size() && bi < 4; ++bi) {
    const DetBranch& br = block.branches[bi];

    // Pool-only branch (Mixed_6a/7a max-pool branch): dispatch
    // max-pool from input directly to br.pool_out.
    if (br.pool_only) {
      if (!max_pool) {
        LOG(ERROR) << "DispatchDetMixedBlock: pool_only branch but max_pool=nil";
        return false;
      }
      MaxPoolDesc mpd = br.max_pool;
      mpd.B = batch_size;
      if (!max_pool->Encode(cb, input_buf, br.pool_out, mpd)) {
        LOG(ERROR) << "max_pool failed (br " << bi << ", " << block.tap_name << ")";
        return false;
      }
      br_outs[bi] = br.pool_out;
      br_c[bi] = br.split_out_C;  // re-used as pool branch C_out
      continue;
    }

    id<MTLBuffer> branch_in = input_buf;
    if (br.has_avg_pool_pre) {
      AvgPoolDesc apd = br.avg_pool;
      apd.B = batch_size;
      if (!avg_pool->Encode(cb, input_buf, br.pool_out, apd)) {
        LOG(ERROR) << "avg_pool failed (br " << bi << ", " << block.tap_name << ")";
        return false;
      }
      branch_in = br.pool_out;
    }

    if (!br.is_split) {
      // Sequential branch.
      for (size_t oi = 0; oi < br.ops.size(); ++oi) {
        const DetBranchOp& op = br.ops[oi];
        ConvDesc cd = op.conv;
        cd.B = batch_size;
        // Folded path: conv writes directly to out_buf (bias + ReLU
        // fused). Unfolded path: conv writes to raw_buf, then BN+ReLU
        // produces out_buf.
        const bool folded = (op.mean == nil);
        id<MTLBuffer> conv_dst = folded ? op.out_buf : op.raw_buf;
        if (!conv_serial->Encode(cb, branch_in, op.w, op.bias,
                                  conv_dst, cd)) {
          LOG(ERROR) << "conv failed (br " << bi << " op " << oi
                     << ", " << block.tap_name << ")";
          return false;
        }
        if (!folded) {
          BnReluDesc bnd{};
          bnd.B = batch_size;
          bnd.H = op.out_H; bnd.W = op.out_W; bnd.C = op.out_C;
          bnd.eps = kBNEpsilon; bnd.relu = true;
          if (!bn_relu->Encode(cb, op.raw_buf, op.mean, op.var, op.beta,
                                op.out_buf, bnd)) {
            LOG(ERROR) << "bn_relu failed (br " << bi << " op " << oi
                       << ", " << block.tap_name << ")";
            return false;
          }
        }
        branch_in = op.out_buf;
      }
      br_outs[bi] = br.ops.back().out_buf;
      br_c[bi] = br.ops.back().out_C;
    } else {
      // Split branch: trunk (sequential ops[0..trunk_size-1]), then 2 parallel
      // ops on trunk_end -> concat into split_concat_out.
      auto encode_op = [&](const DetBranchOp& op, id<MTLBuffer> in) -> bool {
        ConvDesc cd = op.conv;
        cd.B = batch_size;
        const bool folded = (op.mean == nil);
        id<MTLBuffer> conv_dst = folded ? op.out_buf : op.raw_buf;
        if (!conv_serial->Encode(cb, in, op.w, op.bias, conv_dst, cd))
          return false;
        if (!folded) {
          BnReluDesc bnd{};
          bnd.B = batch_size; bnd.H = op.out_H; bnd.W = op.out_W;
          bnd.C = op.out_C; bnd.eps = kBNEpsilon; bnd.relu = true;
          if (!bn_relu->Encode(cb, op.raw_buf, op.mean, op.var, op.beta,
                                op.out_buf, bnd)) return false;
        }
        return true;
      };
      for (int oi = 0; oi < br.trunk_size; ++oi) {
        if (!encode_op(br.ops[oi], branch_in)) return false;
        branch_in = br.ops[oi].out_buf;
      }
      // 2 parallel ops on `branch_in`.
      const DetBranchOp& opa = br.ops[br.trunk_size];
      const DetBranchOp& opb = br.ops[br.trunk_size + 1];
      if (!encode_op(opa, branch_in)) return false;
      if (!encode_op(opb, branch_in)) return false;
      // Intra-branch concat: opa.out_buf + opb.out_buf -> split_concat_out.
      ConcatDesc icd{};
      icd.B = batch_size;
      icd.H = opa.out_H; icd.W = opa.out_W;
      icd.n_branches = 2;
      icd.c_size[0] = opa.out_C;
      icd.c_size[1] = opb.out_C;
      icd.c_size[2] = 0;
      icd.c_size[3] = 0;
      if (!concat->Encode(cb, opa.out_buf, opb.out_buf, nil, nil,
                           br.split_concat_out, icd)) {
        LOG(ERROR) << "intra-branch concat failed (br " << bi << ", "
                   << block.tap_name << ")";
        return false;
      }
      br_outs[bi] = br.split_concat_out;
      br_c[bi] = br.split_out_C;
    }
  }

  // Block-level concat across branches.
  ConcatDesc ccd{};
  ccd.B = batch_size;
  ccd.H = block.H_out; ccd.W = block.W_out;
  ccd.n_branches = static_cast<int>(block.branches.size());
  for (int i = 0; i < 4; ++i) ccd.c_size[i] = (i < ccd.n_branches) ? br_c[i] : 0;
  if (!concat->Encode(cb, br_outs[0], br_outs[1], br_outs[2], br_outs[3],
                       block.concat_out, ccd)) {
    LOG(ERROR) << "block-level concat failed (" << block.tap_name << ")";
    return false;
  }
  return true;
}

}  // namespace deepvariant
