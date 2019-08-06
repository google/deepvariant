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

#include "third_party/nucleus/io/hts_path.h"

#include <stddef.h>

#include "absl/strings/str_cat.h"
#include "third_party/nucleus/platform/types.h"

using absl::StrCat;

namespace nucleus {

const char dflt[] = "";

// Use the default file scheme, unless one is provided.
string fix_path(const std::string &path) {
  size_t i = path.find(':');
  if (i > 0 && i < 20) {
    return path;
  }
  return StrCat(dflt, path);
}

htsFile *hts_open_x(const std::string &fn, const char *mode) {
  string new_path = fix_path(fn);
  return hts_open(new_path.c_str(), mode);
}

htsFile *hts_open_format_x(const std::string &fn, const char *mode,
                           htsFormat *fmt) {
  string new_path = fix_path(fn);
  return hts_open_format(new_path.c_str(), mode, fmt);
}

faidx_t *fai_load3_x(const std::string &fa, const std::string &fai,
                     const std::string &gzi, int flags) {
  string nfa = fix_path(fa);
  string nfai = fix_path(fai);
  string ngzi = fix_path(gzi);
  return fai_load3(nfa.c_str(), nfai.c_str(), ngzi.c_str(), flags);
}

int tbx_index_build_x(const std::string &fn, int min_shift,
                      const tbx_conf_t *conf) {
  string new_path = fix_path(fn);
  return tbx_index_build(new_path.c_str(), min_shift, conf);
}

}  // namespace nucleus
