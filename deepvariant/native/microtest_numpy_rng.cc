// Phase 5.5d/3 microtest — verify NumpyMt19937 + BoundedLemireUint32
// match NumPy 1.24's `np.random.RandomState(seed).randint(...)`
// bit-for-bit on golden vectors captured from
// `google/deepvariant:1.10.0` Docker (numpy 1.24.3, seed 2101079370).
//
// Captured (Docker):
//   randint(0, 1000) ×10:    940, 785, 301,  77, 558, 250, 667, 359, 899, 910
//   randint(0, i+1) i=0..19:   0,   0,   1,   1,   2,   3,   3,   6,   7,   5,
//                               10,   7,   5,   5,   9,   7,   2,   3,   9,   9
//
// Either both blocks PASS or one of the two algorithms (MT or Lemire) is wrong.

#include <cstdio>
#include <cstdint>
#include <vector>

#include "deepvariant/native/numpy_mt19937.h"

int main() {
  using namespace deepvariant::npr;
  int n_fail = 0;

  // Golden vector 1: randint(0, 1000) × 10
  {
    NumpyMt19937 g(2101079370u);
    static const uint32_t expected[10] =
        {940, 785, 301, 77, 558, 250, 667, 359, 899, 910};
    std::printf("Test 1: randint(0, 1000) x10\n");
    bool fail = false;
    for (int i = 0; i < 10; ++i) {
      uint32_t got = RandintU32(g, 1000);
      const char* status = (got == expected[i]) ? "OK" : "FAIL";
      if (got != expected[i]) {
        fail = true;
        std::printf("  [%d] got=%u, expected=%u  %s\n",
                    i, got, expected[i], status);
      }
    }
    if (!fail) std::printf("  → 10/10 match\n");
    else ++n_fail;
  }

  // Golden vector 2: randint(0, i+1) for i = 0..19
  {
    NumpyMt19937 g(2101079370u);
    static const uint32_t expected[20] =
        {0, 0, 1, 1, 2, 3, 3, 6, 7, 5,
         10, 7, 5, 5, 9, 7, 2, 3, 9, 9};
    std::printf("Test 2: randint(0, i+1) for i=0..19\n");
    bool fail = false;
    for (int i = 0; i < 20; ++i) {
      uint32_t got = RandintU32(g, (uint32_t)(i + 1));
      if (got != expected[i]) {
        fail = true;
        std::printf("  [%d] got=%u, expected=%u  FAIL\n", i, got, expected[i]);
      }
    }
    if (!fail) std::printf("  → 20/20 match\n");
    else ++n_fail;
  }

  // Reservoir sample sanity: with k > n, all items kept in order.
  {
    NumpyMt19937 g(2101079370u);
    std::vector<int> items{1, 2, 3, 4, 5};
    auto out = ReservoirSamplePtrs(items, 100, g);
    bool fail = (out.size() != 5);
    for (int i = 0; i < 5 && !fail; ++i) fail |= (*out[i] != items[i]);
    std::printf("Test 3: ReservoirSample k > n preserves order: %s\n",
                fail ? "FAIL" : "OK");
    if (fail) ++n_fail;
  }

  std::printf("\n%d/3 cases FAILED\n", n_fail);
  return n_fail == 0 ? 0 : 1;
}
