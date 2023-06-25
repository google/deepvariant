/*
 * Copyright 2023 Google LLC.
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

// This utility is used to merge phased reads from different shards.
// We can find a consistent phasing if there are reads that overlap multiple
// shards. Please note, that input file must be local, this utility does not
// Google paths.
//
// Usage:
// blaze-bin/learning/genomics/deepvariant/merge_phased_reads_cpp \
// --input_path <Path to sharded tsv file> \
// --output_path <Path to output file> \
// --logtostderr

#include "deepvariant/merge_phased_reads.h"

#include <string>

#include "absl/flags/flag.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "tensorflow/core/platform/logging.h"
#include "re2/re2.h"

ABSL_FLAG(std::string, input_path, "", "Sharded input.");
ABSL_FLAG(std::string, output_path, "", "Output path.");

namespace learning {
namespace genomics {
namespace deepvariant {

absl::StatusOr<ShardedFileSpec>
parse_sharded_file_spec(absl::string_view file_spec) {
  ShardedFileSpec output;
  std::string nshards_str;
  static const absl::string_view kShardSpecPattern =
      R"((.*)\@(\d*[1-9]\d*)(?:\.(.+))?)";
  if (!RE2::FullMatch(file_spec, kShardSpecPattern, &output.basename,
                      &nshards_str, &output.suffix)) {
    // No @<N> found.
    return absl::InvalidArgumentError(
        absl::StrCat("'", file_spec, "' is not a valid sharded file spec."));
  }

  int nshards_int;
  if (!absl::SimpleAtoi(nshards_str, &nshards_int)) {
    return absl::OutOfRangeError(absl::StrCat(
        "'", file_spec, "' is not a valid sharded file spec. '", nshards_str,
        "' is too large to represent as an integer."));
  }
  output.nshards = nshards_int;

  return output;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char** argv) {
  QCHECK(FLAGS_input_path.CurrentValue().empty() ||
      FLAGS_output_path.CurrentValue().empty()) <<
    "ERROR: --input_path and --output_path flags must be set.";

  return 0;
}
