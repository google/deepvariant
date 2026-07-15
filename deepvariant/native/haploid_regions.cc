/*
 * Copyright 2025 Google LLC.
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
 */

#include "deepvariant/native/haploid_regions.h"

#include <zlib.h>

#include <cstdint>
#include <string>
#include <vector>

#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"

namespace deepvariant {

bool LoadParRegions(const std::string& path, ParRegions* out,
                    std::string* err) {
  // gzopen transparently reads both plaintext and gzipped/bgzipped input, so a
  // .bed or .bed.gz both work (matching upstream's htslib BedReader).
  gzFile f = gzopen(path.c_str(), "rb");
  if (f == nullptr) {
    *err = absl::StrCat("cannot open --par_regions_bed: ", path);
    return false;
  }
  // BED lines are short; 64 KiB is far more than any real interval line.
  char buf[1 << 16];
  int lineno = 0;
  bool ok = true;
  while (gzgets(f, buf, sizeof(buf)) != nullptr) {
    ++lineno;
    absl::string_view line(buf);
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
      line.remove_suffix(1);
    }
    if (line.empty() || line[0] == '#') continue;
    if (absl::StartsWith(line, "track") || absl::StartsWith(line, "browser")) {
      continue;
    }
    std::vector<absl::string_view> cols =
        absl::StrSplit(line, absl::ByAnyChar("\t "), absl::SkipEmpty());
    int64_t start, end;
    if (cols.size() < 3 || !absl::SimpleAtoi(cols[1], &start) ||
        !absl::SimpleAtoi(cols[2], &end) || start < 0 || end < start) {
      *err = absl::StrCat("malformed BED line ", lineno, " in ", path, ": ",
                          line);
      ok = false;
      break;
    }
    out->by_contig[std::string(cols[0])].emplace_back(start, end);
  }
  gzclose(f);
  return ok;
}

}  // namespace deepvariant
