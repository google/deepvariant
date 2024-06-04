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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_SSW_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_SSW_H_

#include "src/ssw_cpp.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

// This wrapper library exists solely to provide "string" interface methods,
// since CLIF cannot wrap char*.  Note that we include this wrapper class here
// rather than in SSW, so as to keep things clean for OSS.

class Alignment : public StripedSmithWaterman::Alignment {};

class Filter : public StripedSmithWaterman::Filter {
 public:
  Filter();

  Filter(const bool& pos, const bool& cigar,
         const uint16_t& score, const uint16_t& dis);
};

class Aligner : public StripedSmithWaterman::Aligner {
 public:
  Aligner();

  Aligner(const uint8_t& match_score,
          const uint8_t& mismatch_penalty,
          const uint8_t& gap_opening_penalty,
          const uint8_t& gap_extending_penalty);

  int SetReferenceSequence(const string& reference);

  int Align(const string& query, const Filter& filter, int maskLen,
            Alignment* alignment)
      const;
};


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning



#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_SSW_H_
