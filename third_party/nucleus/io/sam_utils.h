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
 */

#ifndef THIRD_PARTY_NUCLEUS_IO_SAM_UTILS_H_
#define THIRD_PARTY_NUCLEUS_IO_SAM_UTILS_H_

#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/cigar.pb.h"

namespace nucleus {

// This lists the header record type codes defined in section 1.3 of
// https://samtools.github.io/hts-specs/SAMv1.pdf.
//
// Record type tags are found in the header section of SAM files. The
// two-character tag is preceded with @, and indicates the start of a header
// line.
extern const char kSamReadGroupTag[];
extern const char kSamProgramTag[];
extern const char kSamCommentTag[];
extern const char kSamHeaderTag[];
extern const char kSamReferenceSequenceTag[];

// This lists the data field tags defined in section 1.3 of
// https://samtools.github.io/hts-specs/SAMv1.pdf.
//
// Data field tags are found in header lines. Each field tag consists of two
// capital aphabet characters followed by a colon.
//
// These constants are named after k(two-letter-code)Tag, where the
// two-letter-code together with a colon is the contents of the tag strings.
extern const char kCLTag[];
extern const char kCNTag[];
extern const char kDSTag[];
extern const char kDTTag[];
extern const char kFOTag[];
extern const char kGOTag[];
extern const char kIDTag[];
extern const char kKSTag[];
extern const char kLBTag[];
extern const char kLNTag[];
extern const char kPGTag[];
extern const char kPITag[];
extern const char kPLTag[];
extern const char kPMTag[];
extern const char kPNTag[];
extern const char kPPTag[];
extern const char kPUTag[];
extern const char kSMTag[];
extern const char kSNTag[];
extern const char kSOTag[];
extern const char kVNTag[];

// Array mapping CigarUnit_Operation enum to htslib BAM constants.
extern const int kProtoToHtslibCigar[];

// Array mapping htslib BAM constants (in comment) to proto CigarUnit enum
// values.
extern const genomics::v1::CigarUnit_Operation kHtslibCigarToProto[];

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_SAM_UTILS_H_
