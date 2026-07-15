// Phase 5.5 MPSGraph debug walker.
//
// For an all-zeros input, every Inception-v3 stage in the stem produces
// a spatially-constant per-channel output (since each layer's input is
// spatially constant — first layer = relu(bias), and any conv/pool of
// a spatially constant tensor is also spatially constant). After
// Mixed_5b's branches concatenate, the structure stays spatially
// constant for several more stages.
//
// We use that property to localise the first stage where Metal's
// output goes wrong: at every named tap, we sample the channel 0
// value at multiple spatial positions; if they aren't all equal,
// something has injected spatial structure into the all-constant
// input → that's our divergence point.
//
// Also: for stem_s1a we have a closed-form reference and check
// channel-by-channel exactness (32/32 expected on a healthy build).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "deepvariant/native/dv_weights.h"
#include "deepvariant/native/metal_inference.h"

namespace deepvariant {

const DvwTensor* MustGet(const DvwWeights& w, const std::string& name) {
  const auto* t = w.Get(name);
  if (!t) {
    std::fprintf(stderr, "missing tensor: %s\n", name.c_str());
    std::exit(2);
  }
  return t;
}

// Return (B, C, H, W) inferred from the tap's known geometry.
struct TapShape {
  int C, H, W;
};

const std::vector<std::pair<std::string, TapShape>>& TapList() {
  // Shapes are the authoritative shapes produced by the upstream
  // `google/deepvariant:1.10.0` Docker SavedModel forward pass at each
  // named tap, as captured by `tools/conversion/dump_tf_per_layer.py`
  // (see `testdata/reference/per_layer/<tap>.npy`). Metal's MPSGraph
  // builder produces identical shapes (verified 2026-04-28).
  // TapShape is {C, H, W}. With Metal & TF both running NHWC end-to-end
  // (Phase 5.5a fix), per-image tensor sizes are unchanged but the in-
  // memory layout is NHWC. The size check uses C*H*W which is correct
  // either way; downstream compare reads .npy with the matching layout.
  static const std::vector<std::pair<std::string, TapShape>> taps = {
      {"input_nchw", {7, 100, 221}},   // tap kept; NHWC now (size unchanged)
      {"stem_s1a",  {32,  49, 110}},
      {"stem_s2a",  {32,  47, 108}},
      {"stem_s2b",  {64,  47, 108}},
      {"stem_mp3a", {64,  23,  53}},
      {"stem_s3b",  {80,  23,  53}},
      {"stem_s4a", {192,  21,  51}},
      {"stem_mp5a",{192,  10,  25}},
      {"5b",       {256,  10,  25}},
      {"5c",       {288,  10,  25}},
      {"5d",       {288,  10,  25}},
      {"6a",       {768,   4,  12}},
      {"6b",       {768,   4,  12}},
      {"6c",       {768,   4,  12}},
      {"6d",       {768,   4,  12}},
      {"6e",       {768,   4,  12}},
      {"7a",      {1280,   1,   5}},
      {"7b",      {2048,   1,   5}},
      {"7c",      {2048,   1,   5}},
  };
  return taps;
}

// stem_s1a closed-form check for all-zeros input.
// Returns the per-channel-bias-after-ReLU vector for layer 1 (= the
// spatially-constant value of stem_s1a for an all-zero input).
std::vector<float> CheckStemS1a(const DvwWeights& w, MetalInception& inf) {
  constexpr float kEps = 1e-4f;
  const auto* beta = MustGet(w,
      "layer_with_weights-1/beta/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* mean = MustGet(w,
      "layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* var = MustGet(w,
      "layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUE");
  const int O = static_cast<int>(beta->shape[0]);

  std::vector<float> expected(O);
  for (int o = 0; o < O; ++o) {
    const float scale = 1.0f / std::sqrt(var->data[o] + kEps);
    expected[o] = std::max(0.0f, beta->data[o] - mean->data[o] * scale);
  }

  constexpr int B = 1, H = 49, W = 110;
  std::vector<float> input((size_t)B * 100 * 221 * 7, 0.0f);
  std::vector<float> output((size_t)B * O * H * W, 0.0f);
  int per = 0;
  inf.PredictAtTap("stem_s1a", input.data(), B, output.data(), &per);

  int n_match = 0;
  for (int o = 0; o < O; ++o) {
    const size_t idx = (((size_t)0 * O + o) * H + H/2) * W + W/2;
    if (output[idx] == expected[o]) ++n_match;
  }
  std::printf("stem_s1a closed-form: %d/%d channels exact\n", n_match, O);
  return expected;
}

// stem_s2a closed-form check for all-zeros input.
// Input to layer 2 is spatially-constant K_in[c] (the layer-1 bias).
// Layer 2 is conv(3x3 stride 1 valid, in=32, out=32), folded with BN.
// At any *interior* pixel (h,w) of the output, value =
//     relu(b_2[o] + sum_c K_in[c] * sum_{dh,dw} W'_2[o, c, dh, dw])
// where W'_2 is the fold-fused kernel (W * scale_2[o]).
void CheckStemS2a(const DvwWeights& w, MetalInception& inf,
                  const std::vector<float>& k_in) {
  constexpr float kEps = 1e-4f;
  const auto* k = MustGet(w,
      "layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* beta = MustGet(w,
      "layer_with_weights-3/beta/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* mean = MustGet(w,
      "layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* var = MustGet(w,
      "layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUE");
  // kernel shape (3,3,32,32) HWIO; out_dim=32, in_dim=32.
  const int Hk = 3, Wk = 3;
  const int Ik = (int)k->shape[2];
  const int Ok = (int)k->shape[3];

  std::vector<float> expected(Ok);
  for (int o = 0; o < Ok; ++o) {
    const float scale = 1.0f / std::sqrt(var->data[o] + kEps);
    const float bias_o = beta->data[o] - mean->data[o] * scale;
    // Sum over kernel positions and input channels of W * scale * K_in[c].
    float kernel_sum_times_input = 0.0f;
    for (int i = 0; i < Ik; ++i) {
      float kernel_sum_oi = 0.0f;
      for (int h = 0; h < Hk; ++h) {
        for (int wj = 0; wj < Wk; ++wj) {
          const size_t src = ((size_t)h * Wk + wj) * Ik * Ok +
                             (size_t)i * Ok + o;
          kernel_sum_oi += k->data[src];
        }
      }
      kernel_sum_times_input += k_in[i] * (kernel_sum_oi * scale);
    }
    expected[o] = std::max(0.0f, bias_o + kernel_sum_times_input);
  }

  constexpr int B = 1, H = 47, W = 108;
  std::vector<float> input((size_t)B * 100 * 221 * 7, 0.0f);
  std::vector<float> output((size_t)B * Ok * H * W, 0.0f);
  int per = 0;
  inf.PredictAtTap("stem_s2a", input.data(), B, output.data(), &per);

  // Stem_s2a has VALID padding on stride-1 conv → no spatial variation
  // across all positions for spatially-constant input. Sample center
  // pixel for each channel and compare.
  int n_match = 0, n_close = 0;
  float max_diff = 0.0f;
  for (int o = 0; o < Ok; ++o) {
    const size_t idx = (((size_t)0 * Ok + o) * H + H/2) * W + W/2;
    const float metal_v = output[idx];
    const float diff = std::fabs(metal_v - expected[o]);
    max_diff = std::max(max_diff, diff);
    if (metal_v == expected[o]) ++n_match;
    if (diff < 1e-5f) ++n_close;
  }
  std::printf("stem_s2a closed-form: %d/%d exact, %d/%d <1e-5, max diff %.6e\n",
              n_match, Ok, n_close, Ok, max_diff);
}

// At each tap, sample channel 0 at four corners + center. If the
// values differ, the tensor has spatial structure (which it must NOT
// for a uniform all-zeros input). Print the spatial spread.
void WalkTaps(MetalInception& inf) {
  std::printf("tap          C    H    W   ch0[0,0]      ch0[H/2,W/2]   ch0[H-1,W-1]   spread\n");
  std::printf("-----------  ---  ---  ---  ------------  ------------   ------------   -----------\n");
  for (const auto& [name, sh] : TapList()) {
    const int B = 1;
    const size_t total = (size_t)B * sh.C * sh.H * sh.W;
    std::vector<float> input((size_t)B * 100 * 221 * 7, 0.0f);
    std::vector<float> out(total, 0.0f);
    int per = 0;
    if (!inf.PredictAtTap(name, input.data(), B, out.data(), &per)) {
      std::fprintf(stderr, "tap %s failed\n", name.c_str());
      continue;
    }
    if (per != sh.C * sh.H * sh.W) {
      std::printf("%-11s  shape mismatch: per_image=%d, expected C*H*W=%d\n",
                  name.c_str(), per, sh.C * sh.H * sh.W);
      continue;
    }
    auto at = [&](int c, int h, int w) {
      return out[(((size_t)0 * sh.C + c) * sh.H + h) * sh.W + w];
    };
    const float v00 = at(0, 0, 0);
    const float vmid = at(0, sh.H/2, sh.W/2);
    const float vlast = at(0, sh.H-1, sh.W-1);
    // Spread across all spatial positions of channel 0.
    float vmin = v00, vmax = v00;
    for (int h = 0; h < sh.H; ++h) {
      for (int w = 0; w < sh.W; ++w) {
        const float v = at(0, h, w);
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
      }
    }
    std::printf("%-11s  %3d  %3d  %3d  % .6e  % .6e   % .6e   %.3e\n",
                name.c_str(), sh.C, sh.H, sh.W,
                v00, vmid, vlast, vmax - vmin);
  }
}

// ---------------------------------------------------------------------------
// Minimal .npy reader (FP32 little-endian, fortran_order=False).
//
// Format reference: numpy.org/doc/stable/reference/generated/numpy.lib.format
// Header: '\x93NUMPY' + (1B major, 1B minor) + (2B for v1, 4B for v2/3
//   header_len LE) + ASCII header dict (padded with spaces, ends \n) + data.
//
// We only need to extract shape and read the float32 payload.
// ---------------------------------------------------------------------------

struct NpyData {
  std::vector<int> shape;
  std::vector<float> data;
  size_t total = 0;  // = product(shape)
};

bool LoadNpyFp32(const std::string& path, NpyData* out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) {
    std::fprintf(stderr, "npy: cannot open %s\n", path.c_str());
    return false;
  }
  char magic[6];
  f.read(magic, 6);
  if (std::memcmp(magic, "\x93NUMPY", 6) != 0) {
    std::fprintf(stderr, "npy: bad magic in %s\n", path.c_str());
    return false;
  }
  uint8_t major, minor;
  f.read(reinterpret_cast<char*>(&major), 1);
  f.read(reinterpret_cast<char*>(&minor), 1);
  uint32_t header_len;
  if (major == 1) {
    uint16_t hl;
    f.read(reinterpret_cast<char*>(&hl), 2);
    header_len = hl;
  } else {
    uint32_t hl;
    f.read(reinterpret_cast<char*>(&hl), 4);
    header_len = hl;
  }
  std::string header(header_len, '\0');
  f.read(header.data(), header_len);

  // Quick parse: locate "shape" key.
  auto p = header.find("'shape':");
  if (p == std::string::npos) {
    std::fprintf(stderr, "npy: no 'shape' key in header\n");
    return false;
  }
  auto lp = header.find('(', p);
  auto rp = header.find(')', lp);
  if (lp == std::string::npos || rp == std::string::npos) {
    std::fprintf(stderr, "npy: malformed shape\n");
    return false;
  }
  out->shape.clear();
  std::string shape_str = header.substr(lp + 1, rp - lp - 1);
  for (size_t i = 0; i < shape_str.size();) {
    while (i < shape_str.size() &&
           (shape_str[i] == ' ' || shape_str[i] == ',')) {
      ++i;
    }
    if (i >= shape_str.size()) break;
    size_t end = i;
    while (end < shape_str.size() && shape_str[end] >= '0' &&
           shape_str[end] <= '9') {
      ++end;
    }
    if (end == i) break;
    out->shape.push_back(std::stoi(shape_str.substr(i, end - i)));
    i = end;
  }

  // Sanity: descr should be '<f4'.
  if (header.find("'<f4'") == std::string::npos &&
      header.find("'|f4'") == std::string::npos) {
    std::fprintf(stderr, "npy: unsupported dtype in %s (need <f4)\n",
                 path.c_str());
    return false;
  }

  out->total = 1;
  for (int d : out->shape) out->total *= (size_t)d;
  out->data.resize(out->total);
  f.read(reinterpret_cast<char*>(out->data.data()),
         out->total * sizeof(float));
  if (!f) {
    std::fprintf(stderr, "npy: short read on %s\n", path.c_str());
    return false;
  }
  return true;
}

// At each tap, run Metal forward (with the seed-0 input from
// `<ref_dir>/_input.npy`), compare every element to `<ref_dir>/<tap>.npy`,
// print summary stats. Stops at the first tap with max-abs > 1e-3 — that
// is where the structural value bug lives.
int CompareToReference(MetalInception& inf, const std::string& ref_dir) {
  // 1) Load input batch.
  NpyData input;
  if (!LoadNpyFp32(ref_dir + "/_input.npy", &input)) return 1;
  if (input.shape.size() != 4 || input.shape[0] < 1) {
    std::fprintf(stderr, "input shape must be (B, H, W, C); got rank %zu\n",
                 input.shape.size());
    return 1;
  }
  const int B = input.shape[0];
  std::printf("input: shape=(%d, %d, %d, %d), %zu elems\n",
              input.shape[0], input.shape[1], input.shape[2], input.shape[3],
              input.total);
  std::printf("  first 8 vals: ");
  for (int i = 0; i < 8 && i < (int)input.total; ++i) {
    std::printf("%.4f ", input.data[i]);
  }
  std::printf("\n  last 8 vals:  ");
  for (size_t i = input.total > 8 ? input.total - 8 : 0; i < input.total; ++i) {
    std::printf("%.4f ", input.data[i]);
  }
  std::printf("\n");

  // 2) For each tap, run Metal then ULP-diff against ref .npy.
  std::printf("\n%-12s  %-10s  %-12s  %-12s  %-12s  status\n",
              "tap", "n_elems", "max_abs", "mean_abs", "max_rel");
  std::printf("%-12s  %-10s  %-12s  %-12s  %-12s  ------\n",
              "----", "-------", "-------", "--------", "-------");

  int n_ok = 0, n_close = 0, n_diverge = 0;
  for (const auto& [name, sh] : TapList()) {
    NpyData ref;
    const std::string ref_path = ref_dir + "/" + name + ".npy";
    if (!LoadNpyFp32(ref_path, &ref)) continue;
    const size_t total = ref.total;

    std::vector<float> metal_out(total, 0.0f);
    int per = 0;
    if (!inf.PredictAtTap(name, input.data.data(), B, metal_out.data(),
                          &per)) {
      std::printf("%-12s  PredictAtTap failed\n", name.c_str());
      ++n_diverge;
      continue;
    }
    if ((size_t)per * B != total) {
      std::printf("%-12s  size mismatch: per=%d, ref_total=%zu\n",
                  name.c_str(), per, total);
      ++n_diverge;
      continue;
    }

    double max_abs = 0.0, sum_abs = 0.0, max_rel = 0.0;
    for (size_t i = 0; i < total; ++i) {
      const double d = std::fabs((double)metal_out[i] - (double)ref.data[i]);
      sum_abs += d;
      if (d > max_abs) max_abs = d;
      const double denom = std::fabs((double)ref.data[i]);
      if (denom > 1e-6) {
        const double r = d / denom;
        if (r > max_rel) max_rel = r;
      }
    }
    const double mean_abs = sum_abs / (double)total;

    // Threshold rationale: FP32 conv accumulates ~1 ULP / layer
    // (≈ 6e-8 relative) across 188 layers; max-abs at the deepest
    // taps can reach ~5e-3 even when bit-perfect at each step.
    const char* status;
    if (max_abs <= 1e-5) {
      status = "OK";
      ++n_ok;
    } else if (max_abs <= 5e-3) {
      status = "close";
      ++n_close;
    } else {
      status = "DIVERGE";
      ++n_diverge;
    }
    std::printf("%-12s  %-10zu  %-12.6e  %-12.6e  %-12.6e  %s\n",
                name.c_str(), total, max_abs, mean_abs, max_rel, status);
    // Save Metal output for the first divergent tap as .npy for offline
    // Python analysis (FP32 NCHW layout, no header magic — minimal raw
    // dump; counterpart Python loads via np.fromfile).
    if (max_abs > 1e-3 && n_diverge == 1) {
      const std::string raw_path = ref_dir + "/_metal_" + name + ".raw";
      std::ofstream rf(raw_path, std::ios::binary);
      rf.write(reinterpret_cast<const char*>(metal_out.data()),
               total * sizeof(float));
      rf.close();
      std::printf("    raw Metal dump: %s (%zu floats)\n",
                  raw_path.c_str(), total);
    }
    if (max_abs > 1e-3 && n_diverge == 1) {
      // First divergent tap: dump head + tail + stats side-by-side.
      double m_min = 1e30, m_max = -1e30, m_sum = 0;
      double r_min = 1e30, r_max = -1e30, r_sum = 0;
      size_t m_nz = 0, r_nz = 0;
      for (size_t i = 0; i < total; ++i) {
        const float m = metal_out[i], r = ref.data[i];
        if (m < m_min) m_min = m;
        if (m > m_max) m_max = m;
        m_sum += m;
        if (m != 0) ++m_nz;
        if (r < r_min) r_min = r;
        if (r > r_max) r_max = r;
        r_sum += r;
        if (r != 0) ++r_nz;
      }
      std::printf("    Metal: min=%.3f max=%.3f mean=%.3f nonzero=%zu/%zu\n",
                  m_min, m_max, m_sum / total, m_nz, total);
      std::printf("    TF   : min=%.3f max=%.3f mean=%.3f nonzero=%zu/%zu\n",
                  r_min, r_max, r_sum / total, r_nz, total);
      std::printf("    Metal[0..8]:    ");
      for (int i = 0; i < 8 && i < (int)total; ++i) {
        std::printf("%9.3f ", metal_out[i]);
      }
      std::printf("\n    Metal[mid..+8]: ");
      for (int i = 0; i < 8 && (size_t)(total / 2 + i) < total; ++i) {
        std::printf("%9.3f ", metal_out[total / 2 + i]);
      }
      std::printf("\n    Metal[end-8..]: ");
      for (int i = 0; i < 8 && i < (int)total; ++i) {
        std::printf("%9.3f ", metal_out[total - 8 + i]);
      }
      std::printf("\n    TF   [0..8]:    ");
      for (int i = 0; i < 8 && i < (int)total; ++i) {
        std::printf("%9.3f ", ref.data[i]);
      }
      std::printf("\n    TF   [mid..+8]: ");
      for (int i = 0; i < 8 && (size_t)(total / 2 + i) < total; ++i) {
        std::printf("%9.3f ", ref.data[total / 2 + i]);
      }
      std::printf("\n    TF   [end-8..]: ");
      for (int i = 0; i < 8 && i < (int)total; ++i) {
        std::printf("%9.3f ", ref.data[total - 8 + i]);
      }
      std::printf("\n");
    }
  }
  std::printf("\nsummary: %d OK / %d close / %d DIVERGE (of %zu taps)\n",
              n_ok, n_close, n_diverge, TapList().size());
  return n_diverge == 0 ? 0 : 3;
}

int RunDebug(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::fprintf(stderr,
                 "usage:\n"
                 "  %s <wgs.dvw>                 walk + closed-form\n"
                 "  %s <wgs.dvw> <ref_dir>       compare every tap to "
                 "<ref_dir>/<tap>.npy (and use <ref_dir>/_input.npy)\n",
                 argv[0], argv[0]);
    return 2;
  }
  auto w = DvwWeights::Open(argv[1]);
  if (!w) return 1;
  auto inf = MetalInception::Create(argv[1]);
  if (!inf) return 1;

  if (argc == 3) {
    return CompareToReference(*inf, argv[2]);
  }

  auto k_layer1 = CheckStemS1a(*w, *inf);
  CheckStemS2a(*w, *inf, k_layer1);
  std::printf("\n");
  WalkTaps(*inf);
  return 0;
}

}  // namespace deepvariant

int main(int argc, char** argv) {
  return deepvariant::RunDebug(argc, argv);
}
