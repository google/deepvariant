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

#include "third_party/nucleus/io/sam_utils.h"

#include "htslib/sam.h"

namespace nucleus {

const char kSamReadGroupTag[] = "@RG";
const char kSamProgramTag[] = "@PG";
const char kSamCommentTag[] = "@CO";
const char kSamHeaderTag[] = "@HD";
const char kSamReferenceSequenceTag[] = "@SQ";

const char kCLTag[] = "CL:";
const char kCNTag[] = "CN:";
const char kDSTag[] = "DS:";
const char kDTTag[] = "DT:";
const char kFOTag[] = "FO:";
const char kGOTag[] = "GO:";
const char kIDTag[] = "ID:";
const char kKSTag[] = "KS:";
const char kLBTag[] = "LB:";
const char kLNTag[] = "LN:";
const char kPGTag[] = "PG:";
const char kPITag[] = "PI:";
const char kPLTag[] = "PL:";
const char kPMTag[] = "PM:";
const char kPNTag[] = "PN:";
const char kPPTag[] = "PP:";
const char kPUTag[] = "PU:";
const char kSMTag[] = "SM:";
const char kSNTag[] = "SN:";
const char kSOTag[] = "SO:";
const char kVNTag[] = "VN:";

const int kProtoToHtslibCigar[] = {
    // genomics::v1::CigarUnit::OPERATION_UNSPECIFIED,
    BAM_CBACK,
    // genomics::v1::CigarUnit::ALIGNMENT_MATCH,
    BAM_CMATCH,
    // genomics::v1::CigarUnit::INSERT,
    BAM_CINS,
    // genomics::v1::CigarUnit::DELETE,
    BAM_CDEL,
    // genomics::v1::CigarUnit::SKIP,
    BAM_CREF_SKIP,
    // genomics::v1::CigarUnit::CLIP_SOFT,
    BAM_CSOFT_CLIP,
    // genomics::v1::CigarUnit::CLIP_HARD,
    BAM_CHARD_CLIP,
    // genomics::v1::CigarUnit::CPAD,
    BAM_CPAD,
    // genomics::v1::CigarUnit::SEQUENCE_MATCH,
    BAM_CEQUAL,
    // genomics::v1::CigarUnit::SEQUENCE_MISMATCH,
    BAM_CDIFF,
};

// Array mapping htslib BAM constants (in comment) to proto
// genomics::v1::CigarUnit enum values.
const genomics::v1::CigarUnit_Operation kHtslibCigarToProto[] = {
    // #define BAM_CMATCH      0
    genomics::v1::CigarUnit::ALIGNMENT_MATCH,
    // #define BAM_CINS        1
    genomics::v1::CigarUnit::INSERT,
    // #define BAM_CDEL        2
    genomics::v1::CigarUnit::DELETE,
    // #define BAM_CREF_SKIP   3
    genomics::v1::CigarUnit::SKIP,
    // #define BAM_CSOFT_CLIP  4
    genomics::v1::CigarUnit::CLIP_SOFT,
    // #define BAM_CHARD_CLIP  5
    genomics::v1::CigarUnit::CLIP_HARD,
    // #define BAM_CPAD        6
    genomics::v1::CigarUnit::PAD,
    // #define BAM_CEQUAL      7
    genomics::v1::CigarUnit::SEQUENCE_MATCH,
    // #define BAM_CDIFF       8
    genomics::v1::CigarUnit::SEQUENCE_MISMATCH,
    // #define BAM_CBACK       9
    genomics::v1::CigarUnit::OPERATION_UNSPECIFIED,
};

}  // namespace nucleus
