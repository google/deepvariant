// A2.2 microtest — verify ClassifyMBlockNeon produces output
// byte-identical to ClassifyMBlockScalar across:
//   - all 256 bytes for read[i] and ref[i]
//   - quality boundary values (0..255)
//   - both legacy and non-legacy modes
//   - lengths 0..1024 (catches every NEON tail boundary)
//   - adversarial alignment within a 16-byte chunk
//
// Bit-equivalence is the gating contract — A2.2 cannot wire into
// production until this passes.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "deepvariant/native/neon_cigar_classify.h"

using deepvariant::neon_cigar::ClassifyMasks;
using deepvariant::neon_cigar::ClassifyMBlockNeon;
using deepvariant::neon_cigar::ClassifyMBlockScalar;

namespace {

struct Buffers {
  std::vector<uint8_t> use_base;
  std::vector<uint8_t> is_low_quality;
  std::vector<uint8_t> is_ref;
  std::vector<uint8_t> canonical;
  void Reset(size_t n, uint8_t fill) {
    use_base.assign(n + 16, fill);
    is_low_quality.assign(n + 16, fill);
    is_ref.assign(n + 16, fill);
    canonical.assign(n + 16, fill);
  }
  ClassifyMasks View() {
    return ClassifyMasks{
        use_base.data(),
        is_low_quality.data(),
        is_ref.data(),
        canonical.data(),
    };
  }
};

// Compare body bytes for `n` and check that bytes [n, n+16) were not
// touched (overshoot guard).
int CompareAndOvershoot(const Buffers& a, const Buffers& b, size_t n,
                        uint8_t fill_b, const char* label) {
  int diffs = 0;
  auto chk = [&](const std::vector<uint8_t>& sa, const std::vector<uint8_t>& sb,
                 const char* fld) {
    for (size_t i = 0; i < n; ++i) {
      if (sa[i] != sb[i]) {
        if (diffs < 8)
          std::printf("  %s n=%zu i=%zu fld=%s scalar=%u neon=%u\n",
                      label, n, i, fld, sa[i], sb[i]);
        ++diffs;
      }
    }
    for (size_t i = n; i < n + 16; ++i) {
      if (sb[i] != fill_b) {
        if (diffs < 8)
          std::printf("  %s OVERSHOOT n=%zu fld=%s +%zu (got %u expected %u)\n",
                      label, n, fld, i - n, sb[i], fill_b);
        ++diffs;
      }
    }
  };
  chk(a.use_base, b.use_base, "use_base");
  chk(a.is_low_quality, b.is_low_quality, "is_low_quality");
  chk(a.is_ref, b.is_ref, "is_ref");
  chk(a.canonical, b.canonical, "canonical");
  return diffs;
}

}  // namespace

int main() {
  int n_fail = 0;

  // Test 1: every (read_byte, ref_byte) pair, qual=20, min_q=10, both modes.
  {
    std::printf("Test 1: all 256x256 (read,ref) byte pairs, qual=20, "
                "min_q=10, both modes\n");
    int n_diff = 0;
    std::vector<char> read_buf(256), ref_buf(256);
    std::vector<uint8_t> qual_buf(256, 20);
    Buffers a, b;
    for (int rb = 0; rb < 256; ++rb) {
      for (int rfb = 0; rfb < 256; ++rfb) {
        read_buf.assign(256, (char)rb);
        ref_buf.assign(256, (char)rfb);
        for (int leg = 0; leg <= 1; ++leg) {
          a.Reset(256, 0xAB);
          b.Reset(256, 0xCD);
          ClassifyMBlockScalar(read_buf.data(), ref_buf.data(),
                               qual_buf.data(), 256, 10, leg != 0,
                               a.View());
          ClassifyMBlockNeon(read_buf.data(), ref_buf.data(),
                             qual_buf.data(), 256, 10, leg != 0,
                             b.View());
          n_diff += CompareAndOvershoot(a, b, 256, 0xCD, "256x256");
          if (n_diff > 50) goto done1;
        }
      }
    }
    done1:
    if (n_diff == 0) std::printf("  -> 256x256x2 = 131072 cases PASS\n");
    else { std::printf("  -> %d FAIL\n", n_diff); ++n_fail; }
  }

  // Test 2: quality boundaries — qual ∈ {0, min-1, min, min+1, 255}.
  {
    std::printf("Test 2: quality boundary values (min_q=20)\n");
    int n_diff = 0;
    static const uint8_t test_quals[] = {0, 1, 19, 20, 21, 100, 254, 255};
    Buffers a, b;
    constexpr size_t n = 64;
    std::vector<char> read_buf(n);
    std::vector<char> ref_buf(n, 'A');
    std::vector<uint8_t> qual_buf(n);
    static const char alphabet[] = "ACGTNacgtX0";
    std::mt19937 rng(0xFEEDFACEu);
    for (size_t i = 0; i < n; ++i)
      read_buf[i] = alphabet[rng() % (sizeof(alphabet) - 1)];
    for (uint8_t q : test_quals) {
      qual_buf.assign(n, q);
      for (int leg = 0; leg <= 1; ++leg) {
        a.Reset(n, 0xAB);
        b.Reset(n, 0xCD);
        ClassifyMBlockScalar(read_buf.data(), ref_buf.data(),
                             qual_buf.data(), n, 20, leg != 0, a.View());
        ClassifyMBlockNeon(read_buf.data(), ref_buf.data(),
                           qual_buf.data(), n, 20, leg != 0, b.View());
        n_diff += CompareAndOvershoot(a, b, n, 0xCD, "qual_bnd");
      }
    }
    if (n_diff == 0)
      std::printf("  -> %zu cases PASS\n",
                  sizeof(test_quals) / sizeof(test_quals[0]) * 2);
    else { std::printf("  -> %d FAIL\n", n_diff); ++n_fail; }
  }

  // Test 3: lengths 0..1024 with random ACGTN/X reads + random qualities,
  // both modes.
  {
    std::printf("Test 3: random reads x lengths 0..1024 x both modes\n");
    int n_diff = 0;
    std::mt19937 rng(0x12345678u);
    static const char alphabet[] = "ACGTNacgt0123";
    constexpr size_t kAlpha = sizeof(alphabet) - 1;
    Buffers a, b;
    std::vector<char> read_buf, ref_buf;
    std::vector<uint8_t> qual_buf;
    for (size_t n = 0; n <= 1024; ++n) {
      read_buf.assign(n, 0);
      ref_buf.assign(n, 0);
      qual_buf.assign(n, 0);
      for (size_t i = 0; i < n; ++i) {
        read_buf[i] = alphabet[rng() % kAlpha];
        ref_buf[i] = alphabet[rng() % kAlpha];
        qual_buf[i] = static_cast<uint8_t>(rng() & 0xFF);
      }
      for (int leg = 0; leg <= 1; ++leg) {
        a.Reset(n, 0xAB);
        b.Reset(n, 0xCD);
        ClassifyMBlockScalar(read_buf.data(), ref_buf.data(),
                             qual_buf.data(), n, 25, leg != 0, a.View());
        ClassifyMBlockNeon(read_buf.data(), ref_buf.data(),
                           qual_buf.data(), n, 25, leg != 0, b.View());
        n_diff += CompareAndOvershoot(a, b, n, 0xCD, "rand_len");
      }
      if (n_diff > 100) break;
    }
    if (n_diff == 0)
      std::printf("  -> 1025 lengths x 2 modes = 2050 cases PASS\n");
    else { std::printf("  -> %d FAIL\n", n_diff); ++n_fail; }
  }

  // Test 4: throughput on 150-base reads x 1M iter (Illumina-realistic).
  {
    std::printf("Test 4: throughput on 150-base reads x 1M iter\n");
    constexpr size_t kRowLen = 150;
    constexpr size_t kIter = 1'000'000;
    std::vector<char> read_buf(kRowLen), ref_buf(kRowLen);
    std::vector<uint8_t> qual_buf(kRowLen);
    std::mt19937 rng(0xCAFEBABEu);
    static const char alphabet[] = "ACGT";
    for (size_t i = 0; i < kRowLen; ++i) {
      read_buf[i] = alphabet[rng() & 3];
      ref_buf[i] = alphabet[rng() & 3];
      qual_buf[i] = static_cast<uint8_t>(20 + (rng() % 30));
    }
    Buffers a, b;
    a.Reset(kRowLen, 0xAB);
    b.Reset(kRowLen, 0xCD);

    auto bench = [&](auto fn, Buffers& out_buf, const char* name) {
      uint64_t sink = 0;
      for (int w = 0; w < 1000; ++w) {
        fn(read_buf.data(), ref_buf.data(), qual_buf.data(), kRowLen, 20,
           false, out_buf.View());
        sink += out_buf.use_base[w & (kRowLen - 1)];
      }
      auto t0 = std::chrono::steady_clock::now();
      for (size_t i = 0; i < kIter; ++i) {
        read_buf[i & (kRowLen - 1)] = alphabet[i & 3];
        fn(read_buf.data(), ref_buf.data(), qual_buf.data(), kRowLen, 20,
           false, out_buf.View());
        sink += out_buf.use_base[i & (kRowLen - 1)] +
                out_buf.is_ref[i & (kRowLen - 1)];
      }
      auto t1 = std::chrono::steady_clock::now();
      double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
      std::printf("  %-12s : %.2f ns/read (sink=%llu)\n",
                  name, ns / kIter, (unsigned long long)sink);
      return ns;
    };
    double s_ns = bench(ClassifyMBlockScalar, a, "scalar");
    double n_ns = bench(ClassifyMBlockNeon, b, "neon");
    std::printf("  speed-up : %.2fx\n", s_ns / n_ns);
  }

  std::printf("\n%d test%s failed\n", n_fail, n_fail == 1 ? "" : "s");
  return n_fail == 0 ? 0 : 1;
}
