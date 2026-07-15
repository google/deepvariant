// Phase 1 gate: smoke test that nucleus_io static lib compiles and links.
// Instantiates key types without reading actual files.
#include "gtest/gtest.h"
#include "third_party/nucleus/io/hts_verbose.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/gfile.h"

TEST(NucleusIoSmoke, GfileExistsFalse) {
  EXPECT_FALSE(nucleus::Exists("/this/path/does/not/exist"));
}

TEST(NucleusIoSmoke, GlobEmptyOnNonExistentPattern) {
  auto r = nucleus::Glob("/this/does/not/exist/*.bam");
  EXPECT_TRUE(r.empty());
}

TEST(NucleusIoSmoke, HtsVerboseGetSetLevel) {
  // Verify the hts_verbose API compiles and links.
  enum htsLogLevel level = nucleus::HtsGetLogLevel();
  nucleus::HtsSetLogLevel(level);  // round-trip
  SUCCEED();
}
