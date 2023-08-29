/*
 * Copyright 2018 Google LLC.
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

#include "third_party/nucleus/io/reference.h"

#include <stddef.h>
#include <stdlib.h>

#include <algorithm>
#include <utility>

#include "absl/strings/ascii.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "htslib/tbx.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/protos/fasta.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::ReferenceSequence;

// ###########################################################################
//
// GenomeReference code
//
// ###########################################################################

bool GenomeReference::HasContig(const string& contig_name) const {
  const auto& contigs = Contigs();
  return std::any_of(contigs.cbegin(), contigs.cend(),
                     [&](const nucleus::genomics::v1::ContigInfo& contig) {
                       return contig.name() == contig_name;
                     });
}

std::vector<string> GenomeReference::ContigNames() const {
  const auto& contigs = Contigs();
  std::vector<string> keys;
  keys.reserve(contigs.size());
  for (const auto& contig : contigs) {
    keys.push_back(contig.name());
  }
  return keys;
}

StatusOr<const nucleus::genomics::v1::ContigInfo*> GenomeReference::Contig(
    const string& contig_name) const {
  for (const auto& contig : Contigs()) {
    if (contig.name() == contig_name) {
      return &contig;
    }
  }
  return ::nucleus::NotFound(absl::StrCat("Unknown contig ", contig_name));
}

// Note that start and end are 0-based, and end is exclusive. So end
// can go up to the number of bases on contig.
bool GenomeReference::IsValidInterval(const Range& range) const {
  StatusOr<const nucleus::genomics::v1::ContigInfo*> contig_status =
      Contig(range.reference_name());
  if (!contig_status.ok()) return false;
  const int64 n_bases = contig_status.ValueOrDie()->n_bases();
  return range.start() >= 0 && range.start() <= range.end() &&
         range.start() < n_bases && range.end() <= n_bases;
}

// ###########################################################################
//
// IndexedFastaReader code
//
// ###########################################################################

namespace {
// Gets information about the contigs from the fai index faidx.
std::vector<nucleus::genomics::v1::ContigInfo> ExtractContigsFromFai(
    const faidx_t* faidx) {
  int n_contigs = faidx_nseq(faidx);
  std::vector<nucleus::genomics::v1::ContigInfo> contigs(n_contigs);
  for (int i = 0; i < n_contigs; ++i) {
    nucleus::genomics::v1::ContigInfo* contig = &contigs[i];
    const char* name = faidx_iseq(faidx, i);
    CHECK_NE(name, nullptr) << "Name of " << i << " contig in is null";
    contig->set_name(name);
    contig->set_description("");
    contig->set_n_bases(faidx_seq_len(faidx, name));
    CHECK_GE(contig->n_bases(), 0) << "Contig " << name << "Has < 0 bases";
    contig->set_pos_in_fasta(i);
  }
  return contigs;
}
}  // namespace

// Iterable class for traversing all Fasta records in the file.
class IndexedFastaReaderIterable : public GenomeReferenceRecordIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(GenomeReferenceRecord* out) override;

  // Constructor is invoked via IndexedFastaReader::Iterate.
  IndexedFastaReaderIterable(const IndexedFastaReader* reader);
  ~IndexedFastaReaderIterable() override;

 private:
  size_t pos_ = 0;
};

StatusOr<std::unique_ptr<IndexedFastaReader>> IndexedFastaReader::FromFile(
    const string& fasta_path, const string& fai_path, int cache_size_bases) {
  nucleus::genomics::v1::FastaReaderOptions options =
      nucleus::genomics::v1::FastaReaderOptions();
  return FromFile(fasta_path, fai_path, options, cache_size_bases);
}

StatusOr<std::unique_ptr<IndexedFastaReader>> IndexedFastaReader::FromFile(
    const string& fasta_path, const string& fai_path,
    const nucleus::genomics::v1::FastaReaderOptions& options,
    int cache_size_bases) {
  const string gzi = fasta_path + ".gzi";
  faidx_t* faidx = fai_load3_x(fasta_path, fai_path, gzi, 0);
  if (faidx == nullptr) {
    return ::nucleus::NotFound(
        absl::StrCat("could not load fasta and/or fai for fasta ", fasta_path));
  }
  return std::unique_ptr<IndexedFastaReader>(
      new IndexedFastaReader(fasta_path, faidx, options, cache_size_bases));
}

IndexedFastaReader::IndexedFastaReader(
    const string& fasta_path, faidx_t* faidx,
    const nucleus::genomics::v1::FastaReaderOptions& options,
    int cache_size_bases)
    : fasta_path_(fasta_path),
      faidx_(faidx),
      options_(options),
      contigs_(ExtractContigsFromFai(faidx)),
      cache_size_bases_(cache_size_bases),
      small_read_cache_(),
      cached_range_() {}

IndexedFastaReader::~IndexedFastaReader() {
  if (faidx_) {
    NUCLEUS_CHECK_OK(Close());
  }
}

StatusOr<string> IndexedFastaReader::GetBases(const Range& range) const {
  if (faidx_ == nullptr) {
    return ::nucleus::FailedPrecondition(
        "can't read from closed IndexedFastaReader object.");
  }
  if (!IsValidInterval(range))
    return ::nucleus::InvalidArgument(
        absl::StrCat("Invalid interval: ", range.ShortDebugString()));

  if (range.start() == range.end()) {
    // We are requesting an empty string. faidx_fetch_seq does not allow this,
    // so we have to special case it here.
    return string("");
  }

  bool use_cache = (cache_size_bases_ > 0) &&
                   (range.end() - range.start() <= cache_size_bases_);
  Range range_to_fetch;

  if (use_cache) {
    if (cached_range_ && RangeContains(*cached_range_, range)) {
      // Get from cache!
      string result = small_read_cache_.substr(
          range.start() - cached_range_->start(), range.end() - range.start());
      return result;
    } else {
      // Prepare to fetch a sizeable chunk from the FASTA.
      int64 contig_n_bases =
          Contig(range.reference_name()).ValueOrDie()->n_bases();
      range_to_fetch = MakeRange(
          range.reference_name(), range.start(),
          std::min(static_cast<int64>(range.start() + cache_size_bases_),
                   contig_n_bases));
      CHECK(IsValidInterval(range_to_fetch));
    }
  } else {
    range_to_fetch = range;
  }

  // According to htslib docs, faidx_fetch_seq c_name is the contig name,
  // start is the first base (zero-based) to include and end is the last base
  // (zero-based) to include. Len is an output variable returning the length
  // of the fetched region, -2 c_name not present, or -1 for a general error.
  // The returned pointer must be freed. We need to subtract one from our end
  // since end is exclusive in GenomeReference but faidx has an inclusive one.
  int len;
  char* bases =
      faidx_fetch_seq(faidx_, range_to_fetch.reference_name().c_str(),
                      range_to_fetch.start(), range_to_fetch.end() - 1, &len);
  if (len <= 0)
    return ::nucleus::InvalidArgument(
        absl::StrCat("Couldn't fetch bases for ", range.ShortDebugString()));
  string result(bases);
  if (!options_.keep_true_case()) {
    absl::AsciiStrToUpper(&result);
  }
  free(bases);

  if (use_cache) {
    // Update cache.
    small_read_cache_ = result;
    cached_range_ = range_to_fetch;
    // Return the requested substring.
    result = small_read_cache_.substr(0, range.end() - range.start());
  }
  return result;
}

StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>
IndexedFastaReader::Iterate() const {
  return StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>(
      MakeIterable<IndexedFastaReaderIterable>(this));
}

::nucleus::Status IndexedFastaReader::Close() {
  if (faidx_ == nullptr) {
    return ::nucleus::FailedPrecondition("IndexedFastaReader already closed");
  } else {
    fai_destroy(faidx_);
    faidx_ = nullptr;
  }
  return ::nucleus::Status();
}

StatusOr<bool> IndexedFastaReaderIterable::Next(GenomeReferenceRecord* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  const IndexedFastaReader* fasta_reader =
      static_cast<const IndexedFastaReader*>(reader_);
  if (pos_ >= fasta_reader->contigs_.size()) {
    return false;
  }
  const genomics::v1::ContigInfo& contig = fasta_reader->contigs_.at(pos_);
  const string& reference_name = contig.name();
  out->first = reference_name;
  out->second =
      fasta_reader->GetBases(MakeRange(reference_name, 0, contig.n_bases()))
          .ValueOrDie();
  pos_++;
  return true;
}

IndexedFastaReaderIterable::~IndexedFastaReaderIterable() {}

IndexedFastaReaderIterable::IndexedFastaReaderIterable(
    const IndexedFastaReader* reader)
    : Iterable(reader) {}

// ###########################################################################
//
// UnindexedFastaReader code
//
// ###########################################################################

namespace {

// Helper method to get the name in a header line. This function assumes the
// first character is '>'.
absl::string_view GetNameInHeaderLine(absl::string_view line) {
  DCHECK_LT(1, line.size()) << "name must contain more than >";
  size_t space_idx = line.find(' ');
  if (space_idx == string::npos) {
    // No space is found. The name is the entire string after >.
    space_idx = line.size();
  }
  return line.substr(1, space_idx - 1);
}

}  // namespace

// Iterable class for traversing all Fasta records in the file.
class UnindexedFastaReaderIterable : public GenomeReferenceRecordIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(GenomeReferenceRecord* out) override;

  // Constructor is invoked via UnindexedFastaReader::Iterate.
  UnindexedFastaReaderIterable(const UnindexedFastaReader* reader);
  ~UnindexedFastaReaderIterable() override;

 private:
  // If non-empty, contains the name/id in the header line of the next record.
  std::string next_name_;
};

StatusOr<std::unique_ptr<UnindexedFastaReader>> UnindexedFastaReader::FromFile(
    const string& fasta_path) {
  StatusOr<std::unique_ptr<TextReader>> textreader_or =
      TextReader::FromFile(fasta_path);
  NUCLEUS_RETURN_IF_ERROR(textreader_or.status());
  return std::unique_ptr<UnindexedFastaReader>(
      new UnindexedFastaReader(std::move(textreader_or.ValueOrDie())));
}

UnindexedFastaReader::~UnindexedFastaReader() {}

const std::vector<nucleus::genomics::v1::ContigInfo>&
UnindexedFastaReader::Contigs() const {
  LOG(FATAL) << "Unimplemented function invoked : " << __func__;
  return contigs_;
}

StatusOr<string> UnindexedFastaReader::GetBases(const Range& range) const {
  LOG(FATAL) << "Unimplemented function invoked : " << __func__;
  return ::nucleus::Unimplemented("");
}

StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>
UnindexedFastaReader::Iterate() const {
  return StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>(
      MakeIterable<UnindexedFastaReaderIterable>(this));
}

::nucleus::Status UnindexedFastaReader::Close() {
  if (!text_reader_) {
    return ::nucleus::FailedPrecondition("UnindexedFastaReader already closed");
  }
  // Close the file pointer.
  ::nucleus::Status close_status = text_reader_->Close();
  text_reader_ = nullptr;
  return close_status;
}

UnindexedFastaReader::UnindexedFastaReader(
    std::unique_ptr<TextReader> text_reader)
    : text_reader_(std::move(text_reader)) {}

StatusOr<bool> UnindexedFastaReaderIterable::Next(GenomeReferenceRecord* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  DCHECK(out && out->first.empty() && out->second.empty())
      << "out must be default initialized";

  const UnindexedFastaReader* fasta_reader =
      static_cast<const UnindexedFastaReader*>(reader_);
  if (!fasta_reader->text_reader_) {
    return ::nucleus::FailedPrecondition(
        "Cannot iterate a closed UnindexedFastaReader.");
  }
  if (!next_name_.empty()) {
    out->first = next_name_;
    next_name_.clear();
  }
  bool eof = false;
  while (true) {
    // Read one line.
    StatusOr<string> line = fasta_reader->text_reader_->ReadLine();
    if (!line.ok()) {
      if (::nucleus::IsOutOfRange(line.status())) {
        eof = true;
        break;
      }
      return ::nucleus::DataLoss("Failed to parse FASTA");
    }
    std::string l = line.ValueOrDie();

    if (l.empty()) continue;
    // Check if the line is a header or a sequence.
    if (l.at(0) == '>') {
      absl::string_view parsed_name = GetNameInHeaderLine(l);
      if (out->first.empty()) {
        out->first = string(parsed_name);
        continue;
      }
      next_name_ = string(parsed_name);
      return true;
    }
    // Processing a sequence line. If name is absent by now, return an error.
    if (out->first.empty()) {
      return ::nucleus::DataLoss("Name not found in FASTA");
    }
    out->second.append(
        absl::AsciiStrToUpper(absl::StripTrailingAsciiWhitespace(l)));
  }
  if (eof && out->first.empty()) {
    // No more records.
    return false;
  }
  return true;
}

UnindexedFastaReaderIterable::~UnindexedFastaReaderIterable() {}

UnindexedFastaReaderIterable::UnindexedFastaReaderIterable(
    const UnindexedFastaReader* reader)
    : Iterable(reader) {}

// ###########################################################################
//
// InMemoryFastaReader code
//
// ###########################################################################

// Iterable class for traversing all Fasta records in the file.
class FastaFullFileIterable : public GenomeReferenceRecordIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(GenomeReferenceRecord* out) override;

  // Constructor is invoked via InMemoryFastaReader::Iterate.
  FastaFullFileIterable(const InMemoryFastaReader* reader);
  ~FastaFullFileIterable() override;

 private:
  size_t pos_ = 0;
};

// Initializes an InMemoryFastaReader from contigs and seqs.
//
// contigs is a vector describing the "contigs" of this GenomeReference. These
// should include only the contigs present in seqs. A ContigInfo object for a
// contig `chrom` should describe the entire chromosome `chrom` even if the
// corresponding ReferenceSequence only contains a subset of the bases.
//
// seqs is a vector where each element describes a region of the genome we are
// caching in memory and will use to provide bases in the query() operation.
//
// Note that only a single ReferenceSequence for each contig is currently
// supported.
//
// There should be exactly one ContigInfo for each reference_name referred to
// across all ReferenceSequences, and no extra ContigInfos.
StatusOr<std::unique_ptr<InMemoryFastaReader>> InMemoryFastaReader::Create(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<nucleus::genomics::v1::ReferenceSequence>& seqs) {
  std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence> seqs_map;

  for (const auto& seq : seqs) {
    if (seq.region().reference_name().empty() || seq.region().start() < 0 ||
        seq.region().start() > seq.region().end()) {
      return ::nucleus::InvalidArgument(
          absl::StrCat("Malformed region ", seq.region().ShortDebugString()));
    }

    const size_t region_len = seq.region().end() - seq.region().start();
    if (region_len != seq.bases().length()) {
      return ::nucleus::InvalidArgument(
          absl::StrCat("Region size = ", region_len,
                       " not equal to bases.length() ", seq.bases().length()));
    }

    auto insert_result = seqs_map.emplace(seq.region().reference_name(), seq);
    if (!insert_result.second) {
      return ::nucleus::InvalidArgument(absl::StrCat(
          "Each ReferenceSequence must be on a different chromosome but "
          "multiple ones were found on ",
          seq.region().reference_name()));
    }
  }

  return std::unique_ptr<InMemoryFastaReader>(
      new InMemoryFastaReader(contigs, seqs_map));
}

StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>
InMemoryFastaReader::Iterate() const {
  return StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>>(
      MakeIterable<FastaFullFileIterable>(this));
}

StatusOr<string> InMemoryFastaReader::GetBases(const Range& range) const {
  if (!IsValidInterval(range))
    return ::nucleus::InvalidArgument(
        absl::StrCat("Invalid interval: ", range.ShortDebugString()));

  const ReferenceSequence& seq = seqs_.at(range.reference_name());

  if (range.start() < seq.region().start() ||
      range.end() > seq.region().end()) {
    return ::nucleus::InvalidArgument(absl::StrCat(
        "Cannot query range=", range.ShortDebugString(),
        " as this InMemoryFastaReader only has bases in the interval=",
        seq.region().ShortDebugString()));
  }
  const int64 pos = range.start() - seq.region().start();
  const int64 len = range.end() - range.start();
  return seq.bases().substr(pos, len);
}

StatusOr<bool> FastaFullFileIterable::Next(GenomeReferenceRecord* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  const InMemoryFastaReader* fasta_reader =
      static_cast<const InMemoryFastaReader*>(reader_);
  if (pos_ >= fasta_reader->contigs_.size()) {
    return false;
  }
  const string& reference_name = fasta_reader->contigs_.at(pos_).name();
  auto seq_iter = fasta_reader->seqs_.find(reference_name);
  if (seq_iter == fasta_reader->seqs_.end()) {
    return false;
  }
  DCHECK_NE(nullptr, out) << "FASTA record cannot be null";
  out->first = reference_name;
  out->second = seq_iter->second.bases();
  pos_++;
  return true;
}

FastaFullFileIterable::~FastaFullFileIterable() {}

FastaFullFileIterable::FastaFullFileIterable(const InMemoryFastaReader* reader)
    : Iterable(reader) {}

}  // namespace nucleus
