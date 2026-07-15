// Phase 8 / Tier 6.0 — Validate all 11 DetMixedBlocks chained.
//
// Loads stem_mp5a.npy as input, runs the full chain
// 5b → 5c → 5d → 6a → 6b → 6c → 6d → 6e → 7a → 7b → 7c
// using the Det dispatch path, and compares each block's output to
// the corresponding TF reference NPY in /tmp/dv_per_layer/.
//
// Usage:
//   ./microtest_det_inception <wgs.dvw> <ref_dir>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "deepvariant/native/dv_weights.h"
#include "deepvariant/native/metal_avg_pool.h"
#include "deepvariant/native/metal_bn_relu.h"
#include "deepvariant/native/metal_concat.h"
#include "deepvariant/native/metal_conv_serial.h"
#include "deepvariant/native/metal_det_mixed.h"

namespace deepvariant {

struct NpyData {
  std::vector<int> shape;
  std::vector<float> data;
  size_t total = 0;
};

bool LoadNpyFp32(const std::string& path, NpyData* out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  char magic[6];
  f.read(magic, 6);
  if (std::memcmp(magic, "\x93NUMPY", 6) != 0) return false;
  uint8_t major, minor;
  f.read((char*)&major, 1);
  f.read((char*)&minor, 1);
  uint32_t header_len;
  if (major == 1) {
    uint16_t hl;
    f.read((char*)&hl, 2);
    header_len = hl;
  } else {
    uint32_t hl;
    f.read((char*)&hl, 4);
    header_len = hl;
  }
  std::string header(header_len, '\0');
  f.read(header.data(), header_len);
  auto p = header.find("'shape':");
  if (p == std::string::npos) return false;
  auto lp = header.find('(', p);
  auto rp = header.find(')', lp);
  std::string ss = header.substr(lp + 1, rp - lp - 1);
  out->shape.clear();
  for (size_t i = 0; i < ss.size();) {
    while (i < ss.size() && (ss[i] == ' ' || ss[i] == ',')) ++i;
    if (i >= ss.size()) break;
    size_t e = i;
    while (e < ss.size() && ss[e] >= '0' && ss[e] <= '9') ++e;
    if (e == i) break;
    out->shape.push_back(std::stoi(ss.substr(i, e - i)));
    i = e;
  }
  out->total = 1;
  for (int d : out->shape) out->total *= (size_t)d;
  out->data.resize(out->total);
  f.read((char*)out->data.data(), out->total * sizeof(float));
  return (bool)f;
}

void Compare(const char* name, const float* ours, const NpyData& ref) {
  double max_abs = 0.0, sum_abs = 0.0, max_rel = 0.0;
  for (size_t i = 0; i < ref.total; ++i) {
    const double d = std::fabs((double)ours[i] - (double)ref.data[i]);
    sum_abs += d;
    if (d > max_abs) max_abs = d;
    const double denom = std::fabs((double)ref.data[i]);
    if (denom > 1e-6) {
      const double r = d / denom;
      if (r > max_rel) max_rel = r;
    }
  }
  const double mean_abs = sum_abs / (double)ref.total;
  const char* status = (max_abs <= 1e-5) ? "OK"
                     : (max_abs <= 5e-3) ? "close"
                     : "DIVERGE";
  std::printf("%-6s  shape=(%d,%d,%d,%d)  max_abs=%.4e  mean_abs=%.4e  "
              "max_rel=%.4e  %s\n",
              name,
              ref.shape.size() >= 1 ? ref.shape[0] : 0,
              ref.shape.size() >= 2 ? ref.shape[1] : 0,
              ref.shape.size() >= 3 ? ref.shape[2] : 0,
              ref.shape.size() >= 4 ? ref.shape[3] : 0,
              max_abs, mean_abs, max_rel, status);
}

int Run(const std::string& dvw_path, const std::string& ref_dir) {
  // Load TF reference NPYs.
  NpyData input_npy;
  if (!LoadNpyFp32(ref_dir + "/stem_mp5a.npy", &input_npy)) {
    std::fprintf(stderr, "FAIL: cannot load stem_mp5a.npy\n"); return 1;
  }
  const int B = input_npy.shape[0];
  const int H_in = input_npy.shape[1];
  const int W_in = input_npy.shape[2];
  const int C_in = input_npy.shape[3];

  std::printf("Stem output (input): (%d,%d,%d,%d), %zu elems\n",
              B, H_in, W_in, C_in, input_npy.total);

  // Load .dvw + Metal.
  auto dvw_p = DvwWeights::Open(dvw_path);
  if (!dvw_p) { std::fprintf(stderr, "FAIL: open %s\n", dvw_path.c_str()); return 1; }
  const DvwWeights& dvw = *dvw_p;

  id<MTLDevice> device = MTLCreateSystemDefaultDevice();
  if (!device) { std::fprintf(stderr, "FAIL: no Metal device\n"); return 1; }
  id<MTLCommandQueue> queue = [device newCommandQueue];
  auto conv_serial = MetalConvSerial::Create();
  auto bn_relu = MetalBnRelu::Create();
  auto avg_pool = MetalAvgPool::Create();
  auto max_pool = MetalMaxPool::Create();
  auto concat = MetalConcat::Create();
  if (!conv_serial || !bn_relu || !avg_pool || !max_pool || !concat) {
    std::fprintf(stderr, "FAIL: dispatcher\n"); return 1;
  }

  // Build all 11 blocks. Track current geometry through the chain.
  using Builder = bool(*)(id<MTLDevice>, const DvwWeights&, int, int, int, int,
                          DetMixedBlock*);
  struct BlockSpec { Builder fn; const char* name; };
  std::vector<BlockSpec> specs = {
      {BuildDetMixed5b, "5b"},
      {BuildDetMixed5c, "5c"},
      {BuildDetMixed5d, "5d"},
      {BuildDetMixed6a, "6a"},
      {BuildDetMixed6b, "6b"},
      {BuildDetMixed6c, "6c"},
      {BuildDetMixed6d, "6d"},
      {BuildDetMixed6e, "6e"},
      {BuildDetMixed7a, "7a"},
      {BuildDetMixed7b, "7b"},
      {BuildDetMixed7c, "7c"},
  };

  std::vector<DetMixedBlock> blocks(specs.size());
  int H = H_in, W = W_in, C = C_in;
  for (size_t i = 0; i < specs.size(); ++i) {
    if (!specs[i].fn(device, dvw, B, H, W, C, &blocks[i])) {
      std::fprintf(stderr, "FAIL: build %s\n", specs[i].name); return 1;
    }
    H = blocks[i].H_out;
    W = blocks[i].W_out;
    C = blocks[i].C_out;
    std::printf("Built  %s: H=%d W=%d C=%d\n", specs[i].name, H, W, C);
  }

  // Allocate input MTLBuffer + load TF input.
  id<MTLBuffer> input_buf =
      [device newBufferWithBytes:input_npy.data.data()
                          length:input_npy.total * sizeof(float)
                         options:MTLResourceStorageModeShared];

  // Dispatch all 11 blocks chained on a single command buffer.
  id<MTLCommandBuffer> cb = [queue commandBuffer];
  id<MTLBuffer> cur = input_buf;
  for (size_t i = 0; i < blocks.size(); ++i) {
    if (!DispatchDetMixedBlock(cb, conv_serial.get(), bn_relu.get(),
                                avg_pool.get(), max_pool.get(), concat.get(),
                                blocks[i], cur, B)) {
      std::fprintf(stderr, "FAIL: dispatch %s\n", specs[i].name); return 1;
    }
    cur = blocks[i].concat_out;
  }
  [cb commit];
  [cb waitUntilCompleted];

  // Per-block compare to TF reference NPY.
  std::printf("\n%-6s  %-23s  %-13s  %-13s  %-13s  status\n",
              "tap", "shape", "max_abs", "mean_abs", "max_rel");
  for (size_t i = 0; i < blocks.size(); ++i) {
    NpyData ref;
    const std::string ref_path = ref_dir + "/" + specs[i].name + ".npy";
    if (!LoadNpyFp32(ref_path, &ref)) {
      std::fprintf(stderr, "warn: cannot load %s\n", ref_path.c_str());
      continue;
    }
    std::vector<float> our(ref.total, 0.0f);
    std::memcpy(our.data(), [blocks[i].concat_out contents],
                ref.total * sizeof(float));
    Compare(specs[i].name, our.data(), ref);
  }
  return 0;
}

}  // namespace deepvariant

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <wgs.dvw> <ref_dir>\n", argv[0]);
    return 2;
  }
  return deepvariant::Run(argv[1], argv[2]);
}
