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

#include "third_party/nucleus/io/reference_fai.h"

#include <algorithm>

#include "htslib/tbx.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/gtl/optional.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using tensorflow::strings::StrCat;

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

StatusOr<std::unique_ptr<GenomeReferenceFai>> GenomeReferenceFai::FromFile(
    const string& fasta_path, const string& fai_path, int cache_size_bases) {
  const string gzi = fasta_path + ".gzi";
  faidx_t* faidx =
      fai_load3_x(fasta_path.c_str(), fai_path.c_str(), gzi.c_str(), 0);
  if (faidx == nullptr) {
    return tensorflow::errors::NotFound(
        StrCat("could not load fasta and/or fai for fasta ", fasta_path));
  }
  return std::unique_ptr<GenomeReferenceFai>(
      new GenomeReferenceFai(fasta_path, faidx, cache_size_bases));
}

GenomeReferenceFai::GenomeReferenceFai(
    const string& fasta_path, faidx_t* faidx, int cache_size_bases)
    : fasta_path_(fasta_path),
      faidx_(faidx),
      contigs_(ExtractContigsFromFai(faidx)),
      cache_size_bases_(cache_size_bases),
      small_read_cache_(),
      cached_range_() {}

GenomeReferenceFai::~GenomeReferenceFai() {
  if (faidx_) {
    TF_CHECK_OK(Close());
  }
}

StatusOr<string> GenomeReferenceFai::GetBases(const Range& range) const {
  if (faidx_ == nullptr) {
    return tensorflow::errors::FailedPrecondition(
        "can't read from closed GenomeReferenceFai object.");
  }
  if (!IsValidInterval(range))
    return tensorflow::errors::InvalidArgument(
      StrCat("Invalid interval: ", range.ShortDebugString()));

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
            range.start() - cached_range_->start(),
            range.end() - range.start());
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
  char* bases = faidx_fetch_seq(
      faidx_, range_to_fetch.reference_name().c_str(),
      range_to_fetch.start(), range_to_fetch.end() - 1, &len);
  if (len <= 0)
    return tensorflow::errors::InvalidArgument(
        StrCat("Couldn't fetch bases for ", range.ShortDebugString()));
  string result = tensorflow::str_util::Uppercase(bases);
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

string GenomeReferenceFai::Info() const {
  return "GenomeReference backed by htslib FAI index";
}

tensorflow::Status GenomeReferenceFai::Close() {
  if (faidx_ == nullptr) {
    return tensorflow::errors::FailedPrecondition(
        "GenomeReferenceFai already closed");
  } else {
    fai_destroy(faidx_);
    faidx_ = nullptr;
  }
  return tensorflow::Status::OK();
}

}  // namespace nucleus
