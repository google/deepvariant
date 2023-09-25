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

#ifndef THIRD_PARTY_NUCLEUS_IO_MERGE_VARIANTS_H_
#define THIRD_PARTY_NUCLEUS_IO_MERGE_VARIANTS_H_

#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/variant_reader.h"
#include "third_party/nucleus/io/vcf_writer.h"
namespace nucleus {

using nucleus::genomics::v1::Variant;

std::unique_ptr<Variant> CreateRecordFromTemplate(const Variant& t, int start,
                                                  int end,
                                                  const GenomeReference& ref);

// Modifies a variant to include gVCF allele and associated likelihoods.
//
// Modifies the given variant by applying the modification to its alleles and
//    allele-related FORMAT fields.
void TransfromToGvcf(Variant* variant);

// Zero-scales GL to mimic write-then-read.
//
// When writing variants using VcfWriter, GLs are converted to PLs, which is an
// integer format scaled so the most likely genotype has value 0. This function
// modifies the input variant to mimic this transformation of GL -> PL -> GL.
void ZeroScaleGl(Variant* variant);

void MergeAndWriteVariantsAndNonVariants(
    bool only_keep_pass, const std::string& variant_file_path,
    const std::vector<std::string>& non_variant_file_paths,
    const std::string& fasta_path, const std::string& vcf_out_path,
    const std::string& gvcf_out_path,
    const nucleus::genomics::v1::VcfHeader& header,
    bool process_somatic = false);

void MergeAndWriteVariantsAndNonVariants(
    bool only_keep_pass, VariantReader* variant_reader,
    ShardedVariantReader* non_variant_reader, VcfWriter* vcf_writer,
    VcfWriter* gvcf_writer, const GenomeReference& ref,
    bool process_somatic = false);

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_MERGE_VARIANTS_H_
