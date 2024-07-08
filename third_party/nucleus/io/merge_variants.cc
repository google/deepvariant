/*
 * Copyright 2022 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "third_party/nucleus/io/merge_variants.h"

#include <vector>

#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/variant_reader.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

// The alternate allele string for the gVCF "any" alternate allele.
constexpr std::string_view GVCF_ALT_ALLELE = "<*>";

// The genotype likelihood of the gVCF alternate allele for variant calls.
constexpr int _GVCF_ALT_ALLELE_GL = -99;

constexpr std::string_view DEEP_VARIANT_PASS = "PASS";

std::unique_ptr<Variant> CreateRecordFromTemplate(const Variant& t, int start,
                                                  int end,
                                                  const GenomeReference& ref) {
  std::unique_ptr<Variant> v = std::make_unique<Variant>();
  v->MergeFrom(t);
  v->set_start(start);
  v->set_end(end);
  if (start != t.start()) {
    nucleus::genomics::v1::Range range;
    range.set_reference_name(v->reference_name());
    range.set_start(start);
    range.set_end(start + 1);
    *v->mutable_reference_bases() = ref.GetBases(range).ValueOrDie();
  }
  return v;
}

void TransfromToGvcf(Variant* variant) {
  if (std::none_of(
          variant->alternate_bases().cbegin(),
          variant->alternate_bases().cend(),
          [](const std::string& base) { return base == GVCF_ALT_ALLELE; })) {
    variant->mutable_alternate_bases()->Add(GVCF_ALT_ALLELE.data());
    // Add one new GL for het allele/gVCF for each of the other alleles, plus
    // one for the homozygous gVCF allele.
    genomics::v1::VariantCall* call = variant->mutable_calls()->Mutable(0);
    for (size_t i = 0; i < variant->alternate_bases().size() + 1; i++) {
      call->mutable_genotype_likelihood()->Add(_GVCF_ALT_ALLELE_GL);
    }
    genomics::v1::Value iv;
    iv.set_int_value(0);
    if (call->info().contains("AD")) {
      call->mutable_info()->at("AD").mutable_values()->Add(std::move(iv));
    }
    genomics::v1::Value nv;
    nv.set_number_value(0);
    if (call->info().contains("VAF")) {
      call->mutable_info()->at("VAF").mutable_values()->Add(std::move(nv));
    }
  }
}

void ZeroScaleGl(Variant* variant) {
  genomics::v1::VariantCall* call = variant->mutable_calls()->Mutable(0);
  double max_gl = *std::max_element(call->genotype_likelihood().cbegin(),
                                    call->genotype_likelihood().cend());
  for (int i = 0; i < call->genotype_likelihood().size(); ++i) {
    double* gl = call->mutable_genotype_likelihood()->Mutable(i);
    *gl -= max_gl;
  }
}

constexpr int kCacheSize = 300000000;

void MergeAndWriteVariantsAndNonVariants(
    bool only_keep_pass, const std::string& variant_file_path,
    const std::vector<std::string>& non_variant_file_paths,
    const std::string& fasta_path, const std::string& vcf_out_path,
    const std::string& gvcf_out_path,
    const nucleus::genomics::v1::VcfHeader& header,
    const std::vector<nucleus::genomics::v1::Range>& ranges,
    bool process_somatic) {
  // Create VCF and gVCF writers
  nucleus::genomics::v1::VcfWriterOptions writer_options;
  writer_options.set_round_qual_values(true);
  auto writer_or_status =
      nucleus::VcfWriter::ToFile(vcf_out_path, header, writer_options);
  if (!writer_or_status.ok()) {
    LOG(ERROR) << "opening writer failed" << writer_or_status.error_message();
  }
  std::unique_ptr<VcfWriter> vcf_writer =
      std::move(writer_or_status.ValueOrDie());

  writer_or_status =
      nucleus::VcfWriter::ToFile(gvcf_out_path, header, writer_options);
  if (!writer_or_status.ok()) {
    LOG(ERROR) << "opening writer failed" << writer_or_status.error_message();
  }
  std::unique_ptr<VcfWriter> gvcf_writer =
      std::move(writer_or_status.ValueOrDie());

  // Create fasta reader
  std::unique_ptr<IndexedFastaReader> fasta_reader =
      std::move(IndexedFastaReader::FromFile(
                    fasta_path, absl::StrCat(fasta_path, ".fai"), kCacheSize)
                    .ValueOrDie());
  const std::vector<genomics::v1::ContigInfo> contigs = fasta_reader->Contigs();

  absl::flat_hash_map<std::string, uint32_t> contig_index_map;
  for (uint32_t i = 0; i < contigs.size(); i++) {
    contig_index_map[contigs[i].name()] = i;
  }

  // Create reader for variants
  std::unique_ptr<VariantReader> variant_reader =
      VariantReader::Open(variant_file_path, "", contig_index_map);

  // Create reader for non_variants
  std::unique_ptr<ShardedVariantReader> non_variant_reader =
      ShardedVariantReader::Open(non_variant_file_paths, contig_index_map);

  MergeAndWriteVariantsAndNonVariants(only_keep_pass, variant_reader.get(),
                                      non_variant_reader.get(),
                                      vcf_writer.get(), gvcf_writer.get(),
                                      *fasta_reader, ranges, process_somatic);
}

void MergeAndWriteVariantsAndNonVariants(
    bool only_keep_pass, VariantReader* variant_reader,
    ShardedVariantReader* non_variant_reader, VcfWriter* vcf_writer,
    VcfWriter* gvcf_writer, const GenomeReference& ref,
    const std::vector<nucleus::genomics::v1::Range>& ranges,
    bool process_somatic) {
  IndexedVariant variant = variant_reader->GetAndReadNext();

  IndexedVariant nonvariant = non_variant_reader->GetAndReadNext();

  while (variant.variant != nullptr || nonvariant.variant != nullptr) {
    if (!ranges.empty() && nonvariant.variant != nullptr &&
        !nucleus::RangesContainVariant(ranges, *nonvariant.variant)) {
      nonvariant = non_variant_reader->GetAndReadNext();
      continue;
    }
    if (variant.contig_map_index < nonvariant.contig_map_index ||
        (variant.contig_map_index == nonvariant.contig_map_index &&
         variant.variant->end() <= nonvariant.variant->start())) {
      if (!only_keep_pass ||
          (variant.variant->filter().size() == 1 &&
           variant.variant->filter(0) == DEEP_VARIANT_PASS)) {
        if (process_somatic) {
          NUCLEUS_QCHECK_OK(vcf_writer->WriteSomatic(*variant.variant));
        } else {
          NUCLEUS_QCHECK_OK(vcf_writer->Write(*variant.variant));
        }
      }
      ZeroScaleGl(variant.variant.get());
      TransfromToGvcf(variant.variant.get());
      if (process_somatic) {
        NUCLEUS_QCHECK_OK(gvcf_writer->WriteSomatic(*variant.variant));
      } else {
        NUCLEUS_QCHECK_OK(gvcf_writer->Write(*variant.variant));
      }

      variant = variant_reader->GetAndReadNext();
    } else if (nonvariant.contig_map_index < variant.contig_map_index ||
               (nonvariant.contig_map_index == variant.contig_map_index &&
                nonvariant.variant->end() <= variant.variant->start())) {
      if (process_somatic) {
        NUCLEUS_QCHECK_OK(gvcf_writer->WriteSomatic(*nonvariant.variant));
      } else {
        NUCLEUS_QCHECK_OK(gvcf_writer->Write(*nonvariant.variant));
      }

      nonvariant = non_variant_reader->GetAndReadNext();
    } else {
      if (nonvariant.variant->start() < variant.variant->start()) {
        std::unique_ptr<Variant> v = CreateRecordFromTemplate(
            *nonvariant.variant, nonvariant.variant->start(),
            variant.variant->start(), ref);
        if (process_somatic) {
          NUCLEUS_QCHECK_OK(gvcf_writer->WriteSomatic(*v));
        } else {
          NUCLEUS_QCHECK_OK(gvcf_writer->Write(*v));
        }
      }
      if (nonvariant.variant->end() > variant.variant->end()) {
        nonvariant = {.variant = CreateRecordFromTemplate(
                          *nonvariant.variant, variant.variant->end(),
                          nonvariant.variant->end(), ref),
                      .contig_map_index = nonvariant.contig_map_index};
      } else {
        // This non-variant site is subsumed by a Variant. Ignore it.
        nonvariant = non_variant_reader->GetAndReadNext();
      }
    }
  }

  NUCLEUS_QCHECK_OK(vcf_writer->Close());
  NUCLEUS_QCHECK_OK(gvcf_writer->Close());
}

}  // namespace nucleus
