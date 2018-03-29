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

#include "third_party/nucleus/io/reference.h"

#include <algorithm>
#include <numeric>

#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using tensorflow::strings::StrCat;

// ###########################################################################
//
// GenomeReference code
//
// ###########################################################################

int GenomeReference::NContigs() const { return Contigs().size(); }

std::vector<string> GenomeReference::ContigNames() const {
  const auto& contigs = Contigs();
  std::vector<string> keys;
  keys.reserve(contigs.size());
  for (const auto& contig : contigs) {
    keys.push_back(contig.name());
  }
  return keys;
}

bool GenomeReference::HasContig(const string& contig_name) const {
  const auto& contigs = Contigs();
  return std::any_of(contigs.cbegin(), contigs.cend(),
                     [&](const nucleus::genomics::v1::ContigInfo& contig) {
                       return contig.name() == contig_name;
                     });
}

StatusOr<const nucleus::genomics::v1::ContigInfo*> GenomeReference::Contig(
    const string& contig_name) const {
  for (const auto& contig : Contigs()) {
    if (contig.name() == contig_name) {
      return &contig;
    }
  }
  return tensorflow::errors::NotFound(StrCat("Unknown contig ", contig_name));
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

int64 GenomeReference::NTotalBasepairs() const {
  const auto& contigs = Contigs();
  return std::accumulate(
      contigs.cbegin(), contigs.cend(), static_cast<int64>(0),
      [](int64 acc, const nucleus::genomics::v1::ContigInfo& contig) {
        return acc + contig.n_bases();
      });
}

}  // namespace nucleus
