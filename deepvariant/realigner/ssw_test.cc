/*
 * Copyright 2017 Google LLC.
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

#include "deepvariant/realigner/ssw.h"

#include <stdio.h>

#include "absl/strings/str_format.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

// This exercises a bug in gcc 5.4.  internal
// Building libssw with -fno-inline should work around it.
// Updated Oct 2022: Even though we no longer need the workaround, I'm keeping
// this test.
int Gcc54Bug() {
  Aligner a(4, 2, 4, 2);
  Filter f;
  Alignment x;
  a.SetReferenceSequence("tttt");
  int accuracy = a.Align("ttAtt", f, 16, &x);
  if (accuracy != 0) return 1;
  absl::PrintF("accuracy=%d cigar=%s\n", accuracy, x.cigar_string);
  if (x.cigar_string != "2=1I2=") return 1;
  return 0;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main() { return learning::genomics::deepvariant::Gcc54Bug(); }
