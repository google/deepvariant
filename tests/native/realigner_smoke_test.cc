// Phase 1 gate: smoke test that realigner static lib compiles and links.
#include "gtest/gtest.h"
#include "deepvariant/realigner/ssw.h"
#include "deepvariant/realigner/fast_pass_aligner.h"
#include "deepvariant/realigner/window_selector.h"

TEST(RealigerSmoke, SSWAlignmentBasic) {
  // SSW C++ API: set the reference sequence, then align a query.
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  // SetReferenceSequence must be called before Align (4-argument form).
  aligner.SetReferenceSequence("ACGT", 4);
  uint16_t score = aligner.Align("ACGT", filter, &alignment, 15);
  EXPECT_GE(score, 0);
}
