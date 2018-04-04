/*
 * Copyright 2017 Google Inc.
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

// Implementation of allelecounter.h.
#include "deepvariant/allelecounter.h"

#include <algorithm>
#include <cstddef>

#include "deepvariant/utils.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Separator string that will appear between the fragment name and read number
// the string key constructed from a Read with ReadKey().
static constexpr char kFragmentNameReadNumberSeparator[] = "/";

using nucleus::GenomeReference;
using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using tensorflow::StringPiece;
using tensorflow::strings::StrCat;

namespace {

// Resolves UNSPECIFIED allele types to either REFERENCE or SUBSTITITON.
//
// We differentiate between REFERENCE supporting bases and substititions (those
// that are different from the reference base). Rather than inline this logic at
// the point where we create a base Allele from a read, we instead set the type
// to UNSPECIFIED and use this function to determine if the base is REFERENCE
// or a SUBSTITUTION.
AlleleType ResolveAlleleType(const AlleleCount& allele_count,
                             const string& bases, const AlleleType& type) {
  if (type == AlleleType::UNSPECIFIED) {
    DCHECK(bases.size() == 1) << "Expected single base event for UNSPECIFIED at"
                              << allele_count.position().ShortDebugString();
    return allele_count.ref_base() == bases ? AlleleType::REFERENCE
                                            : AlleleType::SUBSTITUTION;
  } else {
    return type;
  }
}

}  // namespace

std::vector<Allele> SumAlleleCounts(const AlleleCount& allele_count) {
  std::map<std::pair<StringPiece, AlleleType>, int> allele_sums;
  for (const auto& entry : allele_count.read_alleles()) {
    ++allele_sums[{entry.second.bases(), entry.second.type()}];
  }
  std::vector<Allele> to_return;
  for (const auto& entry : allele_sums) {
    to_return.push_back(MakeAllele(entry.first.first.ToString(),
                                   entry.first.second, entry.second));
  }

  // Creates a synthetic reference Allele if we saw any reference containing
  // alleles, whose count is tracked (for performance reasons) as an integer
  // in the AlleleCount.ref_supporting_read_count field of the proto. This
  // synthetic allele allows us to provide the same API from this function: a
  // vector of the Alleles observed in allele_count without having to track the
  // read names for reference containing reads, which is very memory-intensive.
  if (allele_count.ref_supporting_read_count() > 0) {
    to_return.push_back(MakeAllele(allele_count.ref_base(),
                                   AlleleType::REFERENCE,
                                   allele_count.ref_supporting_read_count()));
  }

  return to_return;
}

int TotalAlleleCounts(const AlleleCount& allele_count) {
  return allele_count.read_alleles_size() +
         allele_count.ref_supporting_read_count();
}

// Returns true if all the bases in read from offset to offset + len pass
// the quality threshold to be used for generating alleles for our counts.
// offset + len must be less than or equal to the length of the aligned
// sequence of read or a CHECK will fail.
bool CanBasesBeUsed(const nucleus::genomics::v1::Read& read, int offset,
                    int len, const AlleleCounterOptions& options) {
  CHECK_LE(offset + len, read.aligned_quality_size());

  const int min_base_quality = options.read_requirements().min_base_quality();
  for (int i = 0; i < len; i++) {
    if (!nucleus::IsCanonicalBase(read.aligned_sequence()[offset + i])) {
      return false;
    }
    if (read.aligned_quality(offset + i) < min_base_quality) {
      return false;
    }
  }
  return true;
}

AlleleCounter::AlleleCounter(const GenomeReference* const ref,
                             const Range& range,
                             const AlleleCounterOptions& options)
    : ref_(ref), interval_(range), options_(options) {
  // Initialize our counts vector of AlleleCounts with proper position and
  // reference base information. Initially the alleles repeated field is empty.
  const string bases = ref_->GetBases(range).ValueOrDie();
  const int64 len = IntervalLength();
  counts_.reserve(len);
  for (int i = 0; i < len; ++i) {
    AlleleCount allele_count;
    const int64 pos = range.start() + i;
    *(allele_count.mutable_position()) =
        nucleus::MakePosition(range.reference_name(), pos);
    allele_count.set_ref_base(bases.substr(i, 1));
    counts_.push_back(allele_count);
  }
}

string AlleleCounter::RefBases(const int64 rel_start, const int64 len) {
  CHECK_GT(len, 0) << "Length must be >= 1";

  // If our region isn't valid (e.g., it is off the end of the chromosome),
  // return an empty string, otherwise get the actual bases from reference.
  const int abs_start = interval_.start() + rel_start;
  const Range region = nucleus::MakeRange(interval_.reference_name(), abs_start,
                                          abs_start + len);
  if (!ref_->IsValidInterval(region)) {
    return "";
  } else {
    return ref_->GetBases(region).ValueOrDie();
  }
}

string AlleleCounter::GetPrevBase(const Read& read, const int read_offset,
                                  const int interval_offset) {
  CHECK_GE(read_offset, 0) << "read_offset should be 0 or greater";
  if (read_offset == 0) {
    // The read_offset case is here to handle the case where the insertion/
    // deletion/soft_clip is the first cigar element of the read, and there's no
    // previous base in the read, and so we take our previous base from the
    // reference genome instead.
    return RefBases(interval_offset - 1, 1);
  } else {
    // In all other cases we actually take our previous base from the read
    // itself.
    return read.aligned_sequence().substr(read_offset - 1, 1);
  }
}

ReadAllele AlleleCounter::MakeIndelReadAllele(const Read& read,
                                              const int interval_offset,
                                              const int read_offset,
                                              const CigarUnit& cigar) {
  const int op_len = cigar.operation_length();
  const string prev_base = GetPrevBase(read, read_offset, interval_offset);

  if (prev_base.empty() || !nucleus::AreCanonicalBases(prev_base) ||
      (cigar.operation() != CigarUnit::DELETE &&
       !CanBasesBeUsed(read, read_offset, op_len, options_))) {
    // There is no prev_base (we are at the start of the contig), or the bases
    // are unusable, so don't actually add the indel allele.
    return ReadAllele();
  }

  AlleleType type;
  string bases;
  switch (cigar.operation()) {
    case CigarUnit::DELETE:
      type = AlleleType::DELETION;
      bases = RefBases(interval_offset, op_len);
      if (bases.empty()) {
        // We couldn't get the ref bases for the deletion (which can happen if
        // the deletion spans off the end of the contig), so abort now without
        // considering this read any longer. It's rare but such things happen in
        // genomes but they do occur in practice, such as when: (1) the reads
        // spans off the chromosome, but because there's more sequence there
        // (the chromosome isn't complete), which means the read can have
        // whatever CIGAR it likes, which may include a deletion; (2) the
        // chromosome is actually circular, and the aligner is clever enough to
        // know that, and the read's cigar reflect true differences of the read
        // to the alignment at the start of the contig.  Nasty, I know.
        LOG(WARNING) << "Deletion spans off the chromosome for read: "
                     << read.ShortDebugString();
        return ReadAllele();
      }

      if (!nucleus::AreCanonicalBases(bases)) {
        // The reference genome has non-canonical bases that are being deleted.
        // We don't add deletions with non-canonical bases so we return an empty
        // ReadAllele().
        return ReadAllele();
      }

      break;
    case CigarUnit::INSERT:
      type = AlleleType::INSERTION;
      bases = read.aligned_sequence().substr(read_offset, op_len);
      break;
    case CigarUnit::CLIP_SOFT:
      type = AlleleType::SOFT_CLIP;
      bases = read.aligned_sequence().substr(read_offset, op_len);
      break;
    default:
      LOG(FATAL) << "Unexpected cigar operation: " << cigar.DebugString();
  }

  return ReadAllele(interval_offset - 1, StrCat(prev_base, bases), type);
}

void AlleleCounter::AddReadAlleles(const Read& read,
                                   const std::vector<ReadAllele>& to_add) {
  for (size_t i = 0; i < to_add.size(); ++i) {
    const ReadAllele& to_add_i = to_add[i];

    // The read can span beyond and after the interval, so don't add counts
    // outside our interval boundaries.
    if (to_add_i.skip() || to_add_i.position() < 0 ||
        to_add_i.position() >= IntervalLength()) {
      continue;
    }

    // If sequential alleles have the same position, skip the first one. This
    // occurs, for example, when we observe a base at position p on the genome
    // which is enqueued as the ith element of our to_add vector. But the next
    // allele is an indel allele which, because of VCF convention, occurs at
    // position p, is enqueued at i+1 and supercedes the previous base
    // substitution. Resolving these conflicts here allows us to keep the
    // Read => ReadAllele algorithm logic simple.
    if (i + 1 < to_add.size() &&
        to_add_i.position() == to_add[i + 1].position()) {
      continue;
    }

    AlleleCount& allele_count = counts_[to_add_i.position()];
    const AlleleType type =
        ResolveAlleleType(allele_count, to_add_i.bases(), to_add_i.type());

    if (type == AlleleType::REFERENCE) {
      const int prev_count = allele_count.ref_supporting_read_count();
      allele_count.set_ref_supporting_read_count(prev_count + 1);
    } else {
      auto* read_alleles = allele_count.mutable_read_alleles();
      const string key = ReadKey(read);
      const Allele allele = MakeAllele(to_add_i.bases(), type, 1);

      // Naively, there should never be multiple counts for the same read key.
      // We detect such a situation here but only write out a warning. It would
      // be better to have a stronger response (FATAL), but unfortunately we see
      // data in the wild that we need to process that has duplicates.
      if (read_alleles->count(key)) {
        // Not thread safe.
        static int counter = 0;
        if (counter++ < 1) {
          LOG(WARNING) << "Found duplicate read: " << key << " at "
                       << allele_count.position().ShortDebugString();
        }
      }

      (*read_alleles)[key] = allele;
    }
  }
}

void AlleleCounter::Add(const Read& read) {
  // redacted

  const LinearAlignment& aln = read.alignment();
  std::vector<ReadAllele> to_add;
  to_add.reserve(read.aligned_quality_size());
  int read_offset = 0;
  int interval_offset = aln.position().position() - Interval().start();

  for (const auto& cigar_elt : aln.cigar()) {
    const int op_len = cigar_elt.operation_length();
    switch (cigar_elt.operation()) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
        for (int i = 0; i < op_len; ++i) {
          const int base_offset = read_offset + i;
          if (CanBasesBeUsed(read, base_offset, 1, options_)) {
            const string bases = read.aligned_sequence().substr(base_offset, 1);
            to_add.push_back(ReadAllele(interval_offset + i, bases,
                                        AlleleType::UNSPECIFIED));
          }
        }
        read_offset += op_len;
        interval_offset += op_len;
        break;
      case CigarUnit::CLIP_SOFT:
      case CigarUnit::INSERT:
        // Note, by convention VCF insertion/deletion are at the preceding base.
        to_add.push_back(
            MakeIndelReadAllele(read, interval_offset, read_offset, cigar_elt));
        read_offset += op_len;
        // No interval offset change, since an insertion doesn't move us on ref.
        break;
      case CigarUnit::DELETE:
        // By convention VCF insertion/deletion are at the preceding base.
        to_add.push_back(
            MakeIndelReadAllele(read, interval_offset, read_offset, cigar_elt));
        // No read offset change, since a deletion doesn't consume read bases.
        interval_offset += op_len;
        break;
      case CigarUnit::PAD:
      case CigarUnit::SKIP:
        // No read offset change, since a pad/skip don't consume read bases.
        interval_offset += op_len;
        break;
      case CigarUnit::CLIP_HARD:
        break;
      default:
        // Lots of misc. enumerated values from proto that aren't useful such as
        // enumeration values INT_MIN_SENTINEL_DO_NOT_USE_ and
        // OPERATION_UNSPECIFIED.
        break;
    }
  }

  AddReadAlleles(read, to_add);
  ++n_reads_counted_;
}

string AlleleCounter::ReadKey(const Read& read) {
  return StrCat(read.fragment_name(), kFragmentNameReadNumberSeparator,
                read.read_number());
}

void AlleleCounter::Add(const std::vector<Read>& reads) {
  for_each(reads.cbegin(), reads.cend(),
           [this](const Read& read) { this->Add(read); });
}

std::vector<AlleleCountSummary> AlleleCounter::SummaryCounts() const {
  std::vector<AlleleCountSummary> summaries;
  summaries.reserve(counts_.size());
  for (const AlleleCount& allele_count : counts_) {
    AlleleCountSummary summary;
    summary.set_reference_name(allele_count.position().reference_name());
    summary.set_position(allele_count.position().position());
    summary.set_ref_base(allele_count.ref_base());
    summary.set_ref_supporting_read_count(
        allele_count.ref_supporting_read_count());
    summary.set_total_read_count(TotalAlleleCounts(allele_count));
    summary.set_ref_nonconfident_read_count(
        allele_count.ref_nonconfident_read_count());
    summaries.push_back(summary);
  }
  return summaries;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
