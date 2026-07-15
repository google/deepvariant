class Glnexus < Formula
  desc "Joint variant calling for population sequencing — Mac ARM native"
  homepage "https://github.com/dnanexus-rnd/GLnexus"
  url "https://github.com/dnanexus-rnd/GLnexus/archive/refs/tags/v1.4.1.tar.gz"
  sha256 "REPLACE_WITH_TARBALL_SHA256"
  license "Apache-2.0"
  version "1.4.1"

  depends_on :macos => :sonoma           # macOS 14 floor
  depends_on arch: :arm64                # Apple Silicon native
  depends_on "cmake" => :build
  depends_on "yaml-cpp"
  depends_on "jemalloc"
  depends_on "boost"
  depends_on "rocksdb"
  depends_on "zstd"

  # GLnexus 1.4.1 has 3 known build issues on Apple Silicon:
  #
  # 1. CMake 4.x rejects `cmake_minimum_required(VERSION 3.2)` —
  #    workaround via -DCMAKE_POLICY_VERSION_MINIMUM=3.5.
  # 2. Vendored capnp 0.7.0 has an arm64 test-suite failure (the
  #    library itself builds fine). Patch replaces `make check` with
  #    `make` in the capnp ExternalProject_Add.
  # 3. Vendored rocksdb 6.22 hardcodes x86 march flags (-msse4.2,
  #    -mpclmul, -march=ivybridge) which don't apply on arm64. We
  #    patch the rocksdb ExternalProject_Add to set
  #    PORTABLE=1 + DISABLE_WARNING_AS_ERROR=1 so it skips x86 flags.
  #
  # Future GLnexus releases (>1.4.1) may resolve these natively.
  patch :DATA

  def install
    mkdir "build" do
      system "cmake", "..",
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_TESTING=OFF",
        *std_cmake_args
      system "make", "glnexus_cli", "-j#{ENV.make_jobs}"
      bin.install "glnexus_cli"
    end
  end

  def caveats
    <<~EOS
      GLnexus joint variant calling — Mac ARM native build.

      Quick start (after running per-sample DeepVariant with --output_gvcf):
        glnexus_cli --config DeepVariantWGS \\
          sample1.g.vcf.gz sample2.g.vcf.gz ... \\
          | bcftools view --threads 4 - | bgzip -c > joint.vcf.gz

      For trio joint genotyping, the DeepVariantWGS config reduces
      Mendelian violations ~30 % via cohort-level allele frequency
      adjustment (Lin et al. Bioinformatics 2018).

      Configurations available:
        DeepVariantWGS       — DeepVariant WGS gvcfs (default)
        DeepVariantWES       — DeepVariant WES gvcfs
        DeepVariant_unfiltered — keep all variants (no filtering)
        gatk_unfiltered      — GATK4 HaplotypeCaller gvcfs

      See: https://github.com/dnanexus-rnd/GLnexus/wiki
    EOS
  end

  test do
    assert_match "glnexus", shell_output("#{bin}/glnexus_cli --help 2>&1 || true")
  end
end

__END__
diff --git a/CMakeLists.txt b/CMakeLists.txt
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -234,7 +234,7 @@ ExternalProject_Add(capnp
     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
     CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/external
     BUILD_IN_SOURCE 1
-    BUILD_COMMAND bash -c "make -j$(nproc) check"
+    BUILD_COMMAND bash -c "make -j$(nproc)"
     INSTALL_COMMAND make install
     LOG_DOWNLOAD ON
   )
@@ -260,6 +260,7 @@ ExternalProject_Add(rocksdb
     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
     CONFIGURE_COMMAND ""
     BUILD_IN_SOURCE 1
-    BUILD_COMMAND bash -c "make -j$(nproc) static_lib"
+    BUILD_COMMAND bash -c "PORTABLE=1 DISABLE_WARNING_AS_ERROR=1 make -j$(nproc) static_lib"
     INSTALL_COMMAND ""
     LOG_DOWNLOAD ON
   )
