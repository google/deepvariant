// Native realigner — orchestrates upstream's WindowSelector + DeBruijnGraph
// + FastPassAligner to realign reads through assembled haplotypes before
// candidate generation. Without this, ~728 candidate sites that upstream
// finds (positions where reads disagree with the reference and need
// re-alignment) are missing from our pipeline.
//
// Mirrors deepvariant/realigner/realigner.py:Realigner.realign_reads:
//   1. window_selector.select_windows(allele_counter) → list of windows
//   2. for each window: debruijn_graph.build → candidate haplotypes
//   3. assign_reads_to_assembled_regions
//   4. for each assembled_region: fast_pass_aligner.realign_reads
//   5. return realigned reads (preserving unassigned reads as-is)

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/realigner.pb.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace deepvariant {

// Build a RealignerOptions proto with WGS-pipeline defaults (matches
// upstream's flag defaults in realigner.py).
learning::genomics::deepvariant::RealignerOptions DefaultRealignerOptions();

// Realign reads in a region using the assembled-haplotype path. The input
// AlleleCounter must already have the reads added; we only use it to find
// candidate windows. Returns the same number of reads as input — those
// overlapping a candidate window are realigned, the rest are returned
// unchanged.
std::vector<nucleus::genomics::v1::Read> RealignReadsForRegion(
    const std::vector<nucleus::genomics::v1::Read>& reads,
    const nucleus::genomics::v1::Range& region,
    const ::learning::genomics::deepvariant::AlleleCounter& counter,
    const nucleus::GenomeReference& ref_reader,
    const learning::genomics::deepvariant::RealignerOptions& options);

}  // namespace deepvariant
