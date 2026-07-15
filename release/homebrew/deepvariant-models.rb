class DeepvariantModels < Formula
  desc "Models for DeepVariant on Apple Silicon — CoreML, Metal DVW, small-model weights, PON"
  homepage "https://github.com/benjamindemaille/deepvariant"
  version "1.10.0"
  license "BSD-3-Clause"

  # Archive contents (~9.2 GB uncompressed):
  #   *.mlpackage          — CoreML/ANE backend (one per model variant)
  #   *.dvw                — Metal MPSGraph FP32 backend
  #   *_small_weights/     — BNNS-CPU MLP weights (.npy, 6 files each)
  #   deepsomatic_pon/     — Panel-of-Normals VCFs for tumor-only calling:
  #     AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz  (Illumina, ~111 MB)
  #     AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz   (PacBio/ONT, ~254 MB)
  url "https://github.com/benjamindemaille/deepvariant/releases/download/v#{version}/deepvariant-models-#{version}.tar.gz"
  sha256 "REPLACE_WITH_TARBALL_SHA256"

  depends_on :macos => :sonoma
  depends_on arch: :arm64

  def install
    d = share/"deepvariant-models"
    d.install Dir["*.mlpackage"]
    d.install Dir["*.dvw"]
    Dir["*_small_weights"].each { |dir| (d/dir).install Dir["#{dir}/*.npy"] }
    (d/"deepsomatic_pon").install Dir["deepsomatic_pon/*"] if Dir.exist?("deepsomatic_pon")
  end

  def caveats
    <<~EOS
      Models: #{share}/deepvariant-models/

      DeepVariant germline: wgs, wes, pacbio, ont, hybrid, masseq, rnaseq
      DeepTrio:             deeptrio.{wgs,wes,pacbio,ont}_{child,parent}
      DeepSomatic T+N:      deepsomatic.{wgs,wes,pacbio,ont,ffpe_wgs,ffpe_wes}
      DeepSomatic TO:       deepsomatic.*_tumor_only
      Pangenome:            pangenome.wgs

      Each variant ships .mlpackage (ANE/GPU) + .dvw (Metal FP32).
      Small-model weights (*_small_weights/) provided for WGS, PacBio, ONT,
      DeepSomatic WGS, PacBio, ONT, FFPE_WGS variants.
      Panel-of-Normals (auto-selected by model_type — no flag needed):
        Illumina (WGS/WES/FFPE): deepsomatic_pon/AF_ilmn_PON_*.vcf.gz
        PacBio/ONT:              deepsomatic_pon/AF_pacbio_PON_*.vcf.gz

      Override path: export DEEPVARIANT_MODELS_DIR=#{share}/deepvariant-models
    EOS
  end

  test do
    assert_predicate share/"deepvariant-models/wgs.mlpackage", :exist?
    assert_predicate share/"deepvariant-models/wgs.dvw",       :exist?
  end
end
