class Deepvariant < Formula
  desc "Native arm64 macOS DeepVariant — germline/trio/somatic/pangenome + Metal/ANE"
  homepage "https://github.com/benjamindemaille/deepvariant"
  version "1.10.0"
  license "BSD-3-Clause"

  # Bottle-only: arm64 macOS. Build requires htslib/abseil/protobuf/re2/boost
  # + Docker for model conversion — end users get a pre-signed binary.
  bottle do
    root_url "https://github.com/benjamindemaille/deepvariant/releases/download/v#{version}"
    rebuild 0
    sha256 cellar: :any_skip_relocation, arm64_sequoia: "REPLACE_WITH_BOTTLE_SHA256"
    sha256 cellar: :any_skip_relocation, arm64_sonoma:  "REPLACE_WITH_BOTTLE_SHA256"
  end

  depends_on :macos => :sonoma
  depends_on arch: :arm64
  depends_on "htslib"             # bgzip + tabix at runtime
  depends_on "deepvariant-models" # .mlpackage, .dvw, small-model weights, PON

  def install
    bin.install "deepvariant"
    # Multi-call binary symlinks — busybox-style. The deepvariant binary
    # inspects basename(argv[0]) at startup (cli.cc::DetectMultiCall) and
    # dispatches to the right runner. One physical binary, four named
    # entry points; no version-skew risk.
    bin.install_symlink "deepvariant" => "deeptrio"
    bin.install_symlink "deepvariant" => "deepsomatic"
    bin.install_symlink "deepvariant" => "pangenome-aware-deepvariant"
  end

  def caveats
    models = "#{HOMEBREW_PREFIX}/share/deepvariant-models"
    <<~EOS
      Four entry points, one binary (~80 MB). Pick whichever idiom you prefer —
      the canonical `deepvariant <subcommand>` form and the per-tool aliases
      dispatch to the same code:

        deepvariant run                        deepvariant trio
        deepvariant somatic                    deepvariant pangenome
        deeptrio                               deepsomatic
        pangenome-aware-deepvariant

      Quick start (models auto-discovered from the deepvariant-models formula):

        # Germline WGS
        deepvariant run --reads=HG002.bam --ref=GRCh38.fa \\
          --output_vcf=out.vcf --model_type=WGS

        # DeepTrio (or use the canonical: `deepvariant trio ...`)
        deeptrio --reads=child.bam \\
          --reads_parent1=p1.bam --reads_parent2=p2.bam \\
          --ref=GRCh38.fa --model_type=WGS \\
          --output_vcf=child.vcf \\
          --output_vcf_parent1=p1.vcf --output_vcf_parent2=p2.vcf

        # DeepSomatic tumor+normal
        deepsomatic --reads_tumor=tumor.bam --reads_normal=normal.bam \\
          --ref=GRCh38.fa --model_type=WGS --output_vcf=somatic.vcf

        # DeepSomatic tumor-only (with Panel-of-Normals)
        deepsomatic --reads_tumor=tumor.bam --ref=GRCh38.fa \\
          --model_type=WGS_TUMOR_ONLY \\
          --population_vcfs=#{models}/deepsomatic_pon/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz \\
          --output_vcf=tumor_only.vcf

        # Pangenome-aware DV (BAM + GBZ-derived BAM from the upstream
        # Docker preprocessing step — GBZ at runtime is out of scope for v2)
        pangenome-aware-deepvariant --reads=sample.bam \\
          --reads_pangenome=pangenome.bam --ref=GRCh38.fa \\
          --output_vcf=out.vcf

      Get help / version on any entry point:
        deepvariant --version            # version + build SHA
        deepvariant --help               # subcommand list
        deepvariant <subcommand> --help  # flags for that subcommand
        deeptrio --help / --helpfull / --help=<substr>

      ANE acceleration: add --inference_backend=ane_speculate to any command.
      Models directory: #{models}
      Override: export DEEPVARIANT_MODELS_DIR=/custom/path
    EOS
  end

  test do
    # 1. Top-level help is rc=0 and lists all subcommands.
    out = shell_output("#{bin}/deepvariant --help")
    assert_match "Top-level pipelines", out
    assert_match "trio",                out
    assert_match "somatic",             out
    assert_match "pangenome",           out

    # 2. --version reports our tag + upstream version + build SHA.
    ver = shell_output("#{bin}/deepvariant --version")
    assert_match "v2-applesilicon",      ver
    assert_match "DeepVariant #{version}", ver

    # 3. Multi-call symlinks resolve to the same binary and self-identify.
    %w[deeptrio deepsomatic pangenome-aware-deepvariant].each do |alias_name|
      v = shell_output("#{bin}/#{alias_name} --version")
      assert_match alias_name, v
      assert_match "DeepVariant #{version}", v
    end
  end
end
