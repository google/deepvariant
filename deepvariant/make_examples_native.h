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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_MAKE_EXAMPLES_NATIVE_H_
#define LEARNING_GENOMICS_DEEPVARIANT_MAKE_EXAMPLES_NATIVE_H_

#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {


// Generates TensorFlow examples from candidates and reads.
class ExamplesGenerator {
 public:
  explicit ExamplesGenerator(const MakeExamplesOptions& options,
                             bool test_mode = false);

 private:
  friend class ExamplesGeneratorPeer;

  // Generates all pairs of alt alleles and ref.
  std::vector<std::vector<std::string>> AltAlleleCombinations(
      const nucleus::genomics::v1::Variant& variant) const;

  // Creates haplotype by concatenating alt bases and reference from both sides.
  // Haplotype length should be equal pileup image width.
  // ref_start_out and ref_end_out are haplotype's start and end in reference
  // coordinates.
  std::string CreateHaplotype(const nucleus::genomics::v1::Variant& variant,
                              const std::string& alt, int64_t* ref_start_out,
                              int64_t* ref_end_out) const;

  // Make examples config.
  const MakeExamplesOptions options_;

  std::unique_ptr<nucleus::GenomeReference> ref_reader_;

  // Half width of pileup image.
  int half_width_;
};

// Helper class to allow unit testing of some private methods.
class ExamplesGeneratorPeer {
 public:
  // Calls private AltAlleleCombinations.
  static std::vector<std::vector<std::string>> CallAltAlleleCombinations(
      const ExamplesGenerator& generator,
      const nucleus::genomics::v1::Variant& variant) {
    return generator.AltAlleleCombinations(variant);
  }

  static std::string CallCreateHaplotype(
      const ExamplesGenerator& generator,
      const nucleus::genomics::v1::Variant& variant, const std::string& alt,
      int64_t* ref_start_out, int64_t* ref_end_out) {
    return generator.CreateHaplotype(variant, alt, ref_start_out, ref_end_out);
  }

  static void SetRefReader(
      ExamplesGenerator& generator,
      std::unique_ptr<nucleus::GenomeReference> ref_reader) {
    generator.ref_reader_ = std::move(ref_reader);
  }
};

nucleus::genomics::v1::Range MakeRange(const std::string& ref_name,
                                       int64_t start, int64_t end);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_MAKE_EXAMPLES_NATIVE_H_
