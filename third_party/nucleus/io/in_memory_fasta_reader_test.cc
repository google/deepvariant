/*
 * Copyright 2018 Google Inc.
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

#include "third_party/nucleus/io/in_memory_fasta_reader.h"

#include <string>

#include "third_party/nucleus/util/utils.h"

#include "tensorflow/core/platform/test.h"

namespace nucleus {

namespace {

// Helper method to create a test sequence.
void CreateTestSeq(std::vector<genomics::v1::ContigInfo>* contigs,
                   std::vector<genomics::v1::ReferenceSequence>* seqs,
                   const string& name, const int pos_in_fasta,
                   const int range_start, const int range_end,
                   const string& bases) {
  DCHECK(pos_in_fasta >= 0 && pos_in_fasta < contigs->size());
  genomics::v1::ContigInfo* contig = &contigs->at(pos_in_fasta);
  contig->set_name(name);
  contig->set_pos_in_fasta(pos_in_fasta);
  contig->set_n_bases(range_end - range_start);
  genomics::v1::ReferenceSequence* seq = &seqs->at(pos_in_fasta);
  seq->mutable_region()->set_reference_name(name);
  seq->mutable_region()->set_start(range_start);
  seq->mutable_region()->set_end(range_end);
  seq->set_bases(bases);
}

}  // namespace

TEST(InMemoryFastaReaderTest, TestIterate) {
  int kNum = 3;
  std::vector<genomics::v1::ContigInfo> contigs(kNum);
  std::vector<genomics::v1::ReferenceSequence> seqs(kNum);
  CreateTestSeq(&contigs, &seqs, "Chr1", 0, 0, 1, "A");
  CreateTestSeq(&contigs, &seqs, "Chr2", 1, 4, 6, "CG");
  CreateTestSeq(&contigs, &seqs, "Chr3", 2, 10, 15, "AATTC");

  std::unique_ptr<InMemoryFastaReader> reader =
      std::move(InMemoryFastaReader::Create(contigs, seqs).ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  StatusOr<bool> status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr1", r.first);
  EXPECT_EQ("A", r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr2", r.first);
  EXPECT_EQ("CG", r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr3", r.first);
  EXPECT_EQ("AATTC", r.second);

  // Reading beyond the file fails.
  status = iterator->Next(&r);
  EXPECT_FALSE(status.ValueOrDie());
}

}  // namespace nucleus
