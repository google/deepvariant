// Phase 8 / Tier 6.0 — Deterministic Inception block dispatch.
//
// Each of the 11 Mixed_X blocks (5b, 5c, 5d, 6a, 6b-6e, 7a, 7b-7c) is
// encoded as a DetMixedBlock with its branches' raw conv weights +
// BN params + intermediate MTLBuffers. Dispatch fans out branches in
// parallel onto a single MTLCommandBuffer, then concats along the
// channel axis.
//
// Bypasses MPSGraph entirely → output is bit-deterministic across
// reduction orders (per-thread sequential FMA via MetalConvSerial).

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "deepvariant/native/dv_weights.h"
#include "deepvariant/native/metal_avg_pool.h"
#include "deepvariant/native/metal_bn_relu.h"
#include "deepvariant/native/metal_concat.h"
#include "deepvariant/native/metal_conv_serial.h"

#ifdef __OBJC__
@protocol MTLDevice;
@protocol MTLCommandBuffer;
@protocol MTLBuffer;
#endif

namespace deepvariant {

struct DetBranchOp {
#ifdef __OBJC__
  ConvDesc conv;
  id<MTLBuffer> w;          // raw HWIO kernel
  id<MTLBuffer> bias;       // all-zero (unfolded BN)
  id<MTLBuffer> mean;
  id<MTLBuffer> var;
  id<MTLBuffer> beta;
  id<MTLBuffer> raw_buf;    // post-conv pre-BN
  id<MTLBuffer> out_buf;    // post-BN+ReLU
#endif
  int out_H = 0, out_W = 0, out_C = 0;
};

// One branch within a Mixed block.
//
// Sequential branch (`is_split == false`):
//   ops[0] -> ops[1] -> ... -> ops[N-1]
//   Branch output buffer = ops.back().out_buf
//
// Split branch (`is_split == true`):
//   ops[0..trunk_size-1] form a sequential trunk.
//   ops[trunk_size] and ops[trunk_size+1] both consume the trunk's
//   final output and run in parallel (e.g. 1×3 and 3×1 in Mixed_7b/7c).
//   Branch output = concat(ops[trunk_size].out_buf, ops[trunk_size+1].out_buf)
//   stored in split_concat_out.
struct DetBranch {
  bool has_avg_pool_pre = false;
  bool has_max_pool_pre = false;
  AvgPoolDesc avg_pool{};
  MaxPoolDesc max_pool{};
#ifdef __OBJC__
  id<MTLBuffer> pool_out;          // input to first op when has_*_pool_pre
#endif
  bool pool_only = false;          // branch = pool only, no convs (Mixed_6a/7a max-pool branch)

  std::vector<DetBranchOp> ops;

  // Split-branch fields (Mixed_7b/7c only):
  bool is_split = false;
  int trunk_size = 0;              // # ops in sequential trunk before split
#ifdef __OBJC__
  id<MTLBuffer> split_concat_out;
#endif
  int split_out_C = 0;             // c_size for block-level concat
};

struct DetMixedBlock {
  std::string tap_name;
  int B = 0;
  int H_in = 0, W_in = 0, C_in = 0;
  int H_out = 0, W_out = 0, C_out = 0;
  std::vector<DetBranch> branches;
#ifdef __OBJC__
  id<MTLBuffer> concat_out;
#endif
};

// Builders for each block type. Each loads the .dvw weights and
// allocates per-branch intermediate buffers + output buffer.
//
// Returns false on weight-load or alloc failure.
#ifdef __OBJC__
bool BuildDetMixed5b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed5c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed5d(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed6a(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed6b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed6c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed6d(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed6e(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed7a(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed7b(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);
bool BuildDetMixed7c(id<MTLDevice> device, const DvwWeights& dvw,
                     int max_B, int H_in, int W_in, int C_in,
                     DetMixedBlock* out);

// Dispatch one block onto `cb`. Reads from input_buf, writes
// concatenated output to block.concat_out. `max_pool` may be null if
// no block in the chain uses a max-pool branch (only Mixed_6a/7a do).
bool DispatchDetMixedBlock(id<MTLCommandBuffer> cb,
                           MetalConvSerial* conv_serial,
                           MetalBnRelu* bn_relu,
                           MetalAvgPool* avg_pool,
                           MetalMaxPool* max_pool,
                           MetalConcat* concat,
                           const DetMixedBlock& block,
                           id<MTLBuffer> input_buf,
                           int batch_size);
#endif

}  // namespace deepvariant
