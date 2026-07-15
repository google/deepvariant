// PoC for option-2 borderline-only-CPU rerun. Computes stem_s1a
// (first conv + BN + ReLU of Inception-v3) on CPU via strict-scalar
// FP32, single-thread, sequential reduction; compares against the
// TF reference dumped from `dump_tf_per_layer.py` inside the
// google/deepvariant:1.10.0 Docker.
//
// Outcome decides whether to invest in a full BNNS-CPU Inception-v3
// for borderline-only re-evaluation:
//
//   - 0 ULP / element  → full bit-exact path achievable (continue)
//   - ≤ 2 ULP / element → close enough; per-layer drift bounded
//                          by sqrt(188) × 2 ≈ 27 ULP at softmax,
//                          which is FILTER-robust; continue
//   - ≫ 2 ULP / element → scalar arm64 ≠ TF AVX-512 at the layer
//                          level; need scalar Docker reference
//                          re-capture before continuing (master
//                          plan §5.5g.0c).

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

#include "deepvariant/native/dv_weights.h"

namespace {

constexpr float kBNEpsilon = 1e-3f;  // Keras default (NOT 1e-4)

bool ReadRawFloats(const std::string& path, std::vector<float>* out,
                   size_t expected) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f) {
    std::fprintf(stderr, "  open %s failed\n", path.c_str());
    return false;
  }
  const std::streamsize sz = f.tellg();
  if ((size_t)sz != expected * sizeof(float)) {
    std::fprintf(stderr, "  %s: %lld bytes, expected %zu floats (%zu bytes)\n",
                 path.c_str(), (long long)sz, expected, expected * sizeof(float));
    return false;
  }
  f.seekg(0);
  out->resize(expected);
  f.read(reinterpret_cast<char*>(out->data()), sz);
  return f.good();
}

uint32_t Ulp(float a, float b) {
  if (a == b) return 0;
  if (std::isnan(a) || std::isnan(b)) return UINT32_MAX;
  uint32_t ai, bi;
  std::memcpy(&ai, &a, 4);
  std::memcpy(&bi, &b, 4);
  // Handle sign asymmetry: use sign-magnitude → 2's-complement-like
  if (ai & 0x80000000u) ai = 0x80000000u - (ai & 0x7FFFFFFFu);
  if (bi & 0x80000000u) bi = 0x80000000u - (bi & 0x7FFFFFFFu);
  return (ai > bi) ? (ai - bi) : (bi - ai);
}

}  // namespace

int main(int argc, char** argv) {
  using namespace deepvariant;

  const std::string dvw_path = (argc > 1)
      ? argv[1]
      : "/Users/benjamin/deepvariant/validation/work/wgs.dvw";
  const std::string ref_dir = "/tmp/dv_per_layer";

  std::printf("=== A/B test: scalar BNNS-CPU stem_s1a vs TF Docker reference ===\n");
  std::printf("dvw         : %s\n", dvw_path.c_str());
  std::printf("ref_dir     : %s\n", ref_dir.c_str());

  // 1) Open weights bundle.
  auto W = DvwWeights::Open(dvw_path);
  if (!W) {
    std::fprintf(stderr, "FATAL: cannot open .dvw at %s\n", dvw_path.c_str());
    return 2;
  }

  // 2) Pull layer-0 (conv2d) and layer-1 (BN) tensors.
  const auto* k = W->Get(
      "layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* beta = W->Get(
      "layer_with_weights-1/beta/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* mean_t = W->Get(
      "layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUE");
  const auto* var_t = W->Get(
      "layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUE");
  if (!k || !beta || !mean_t || !var_t) {
    std::fprintf(stderr, "FATAL: missing weight tensor in .dvw\n");
    return 2;
  }
  if (k->shape.size() != 4 || k->shape[0] != 3 || k->shape[1] != 3 ||
      k->shape[2] != 7 || k->shape[3] != 32) {
    std::fprintf(stderr, "FATAL: unexpected kernel shape\n");
    return 2;
  }
  std::printf("kernel HWIO : %u,%u,%u,%u\n",
              k->shape[0], k->shape[1], k->shape[2], k->shape[3]);

  // 3) Load TF reference input + stem_s1a.
  std::vector<float> input;
  if (!ReadRawFloats(ref_dir + "/_input.raw", &input, 1u * 100 * 221 * 7))
    return 2;
  std::vector<float> ref;
  if (!ReadRawFloats(ref_dir + "/stem_s1a.raw", &ref, 1u * 49 * 110 * 32))
    return 2;
  std::printf("input       : (1,100,221,7) loaded ; ref (1,49,110,32) loaded\n");

  // 4) Strict-scalar FP32 conv (NHWC input, HWIO kernel, stride 2 valid)
  //    → BN → ReLU. Sequential reduction order:
  //    for o in 0..32:
  //      acc = 0
  //      for kh in 0..3:
  //        for kw in 0..3:
  //          for i in 0..7:
  //            acc += input[h*2+kh, w*2+kw, i] * kernel[kh, kw, i, o]
  //      acc = (acc - mean[o]) * (1/sqrt(var[o] + eps)) + beta[o]
  //      acc = max(0, acc)
  //
  // No SIMD, no FMA, no parallel reduction. The C++ compiler under
  // -fno-fast-math (default in CMakeLists.txt) emits sequential
  // mul + add operations matching IEEE 754 strictly.
  constexpr int H_in = 100, W_in = 221, C_in = 7;
  constexpr int H_out = 49, W_out = 110, C_out = 32;
  constexpr int Kh = 3, Kw = 3, Sh = 2, Sw = 2;

  std::vector<float> scale(C_out), offset(C_out);
  for (int o = 0; o < C_out; ++o) {
    scale[o] = 1.0f / std::sqrt(var_t->data[o] + kBNEpsilon);
    offset[o] = beta->data[o] - mean_t->data[o] * scale[o];
  }

  std::vector<float> out(H_out * W_out * C_out);
  for (int h = 0; h < H_out; ++h) {
    for (int w = 0; w < W_out; ++w) {
      for (int o = 0; o < C_out; ++o) {
        float acc = 0.0f;
        for (int kh = 0; kh < Kh; ++kh) {
          const int ih = h * Sh + kh;
          if (ih >= H_in) continue;  // valid padding: skip OOB
          for (int kw = 0; kw < Kw; ++kw) {
            const int iw = w * Sw + kw;
            if (iw >= W_in) continue;
            for (int i = 0; i < C_in; ++i) {
              const float x = input[(ih * W_in + iw) * C_in + i];
              const float wt = k->data[((kh * Kw + kw) * C_in + i) * C_out + o];
              acc += x * wt;
            }
          }
        }
        // BN with raw conv output (separate path, mirrors TF's
        // unfolded Conv→BN→ReLU graph in the frozen reference)
        const float bn = acc * scale[o] + offset[o];
        out[(h * W_out + w) * C_out + o] = bn > 0 ? bn : 0;
      }
    }
  }

  // 5) Compare to TF reference element-by-element.
  uint32_t max_ulp = 0, sum_ulp = 0;
  uint32_t worst_idx = 0;
  double max_abs = 0.0, sum_abs = 0.0;
  size_t n_zero_match = 0, n_zero_diff = 0;
  size_t n = out.size();
  for (size_t i = 0; i < n; ++i) {
    const uint32_t u = Ulp(out[i], ref[i]);
    if (u > max_ulp) { max_ulp = u; worst_idx = i; }
    sum_ulp += (u > 1u << 24 ? 1u << 24 : u);  // saturate for sum
    const double abs_d = std::fabs((double)out[i] - (double)ref[i]);
    if (abs_d > max_abs) max_abs = abs_d;
    sum_abs += abs_d;
    if (out[i] == 0.0f && ref[i] == 0.0f) ++n_zero_match;
    if ((out[i] == 0.0f) != (ref[i] == 0.0f)) ++n_zero_diff;
  }
  std::printf("\n=== Element-wise diff ours-vs-TF on %zu output elements ===\n", n);
  std::printf("  max ULP             : %u  (at idx %u: ours=%.7g, ref=%.7g)\n",
              max_ulp, worst_idx, out[worst_idx], ref[worst_idx]);
  std::printf("  mean ULP            : %.3f\n", (double)sum_ulp / n);
  std::printf("  max abs diff        : %.4g\n", max_abs);
  std::printf("  mean abs diff       : %.4g\n", sum_abs / n);
  std::printf("  zero-on-both        : %zu (%.1f %%)\n",
              n_zero_match, 100.0 * n_zero_match / n);
  std::printf("  zero-mismatch (ReLU bound): %zu\n", n_zero_diff);

  // 6) Verdict.
  std::printf("\n=== Verdict ===\n");
  if (max_ulp == 0) {
    std::printf("  EXCELLENT: scalar BNNS-CPU is BIT-EXACT vs TF reference.\n"
                "  → Full BNNS-CPU big-CNN path is feasible; continue.\n");
    return 0;
  }
  if (max_ulp <= 2) {
    std::printf("  GOOD: scalar BNNS-CPU within %u ULP of TF reference.\n"
                "  → 188-layer cumulative drift bound: sqrt(188)*%u = ~%u ULP.\n"
                "  → FILTER-robust at GQ=20 boundary; continue full path.\n",
                max_ulp, max_ulp, (unsigned)(13 * max_ulp));
    return 0;
  }
  std::printf("  CAUTION: scalar BNNS-CPU drifts %u ULP from TF reference.\n"
              "  → arm64 scalar order != x86 oneDNN AVX-512 reduction tree.\n"
              "  → Re-capture Docker reference in scalar mode\n"
              "    (TF_DISABLE_MKL=1, TF_NUM_INTRAOP_THREADS=1)\n"
              "    before investing in full BNNS-CPU forward pass.\n",
              max_ulp);
  return 1;
}
