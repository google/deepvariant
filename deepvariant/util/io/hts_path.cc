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

#include "third_party/nucleus/io/hts_path.h"
#include <string>
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "tensorflow/core/lib/strings/strcat.h"

using tensorflow::strings::StrCat;
using tensorflow::string;

namespace nucleus {

const char dflt[] = "";

// Use the default file scheme, unless one is provided.
string fix_path(const char *path) {
  string s(path);
  size_t i = s.find(':');
  if (i > 0 && i < 20) {
    return s;
  }
  return StrCat(dflt, s);
}

htsFile *hts_open_x(const char *path, const char *mode) {
  string new_path = fix_path(path);
  return hts_open(new_path.c_str(), mode);
}

faidx_t *fai_load3_x(const char *fa, const char *fai, const char *gzi,
                     int flags) {
  string nfa = fix_path(fa);
  string nfai = fix_path(fai);
  string ngzi = fix_path(gzi);
  return fai_load3(fa ? nfa.c_str() : nullptr, fai ? nfai.c_str() : nullptr,
                   gzi ? ngzi.c_str() : nullptr, flags);
}

}  // namespace nucleus
