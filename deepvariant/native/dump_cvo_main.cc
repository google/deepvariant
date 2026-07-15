// Tiny dumper: read a TFRecord stream of CallVariantsOutput protos and
// print one line per record:
//     <chrom>\t<start1>\t<ref>\t<alt0>\t<alt1>...\t<argmax>
// Used during realigner/AlleleCounter parity work to diff our candidate
// set against upstream's. Not shipped in releases.
//
// Usage: dump_cvo <tfrecord-path>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

#include "deepvariant/native/tfrecord.h"
#include "deepvariant/protos/deepvariant.pb.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <cvo.tfrecord>\n", argv[0]);
    return 2;
  }
  auto rdr = deepvariant::TFRecordReader::New(argv[1]);
  if (!rdr) {
    std::fprintf(stderr, "failed to open %s\n", argv[1]);
    return 1;
  }
  long total = 0;
  while (rdr->GetNext()) {
    learning::genomics::deepvariant::CallVariantsOutput cvo;
    if (!cvo.ParseFromString(rdr->record())) continue;
    const auto& v = cvo.variant();
    int argmax = 0;
    double best = -1.0;
    for (int i = 0; i < cvo.genotype_probabilities_size(); ++i) {
      if (cvo.genotype_probabilities(i) > best) {
        best = cvo.genotype_probabilities(i);
        argmax = i;
      }
    }
    std::cout << v.reference_name() << '\t' << (v.start() + 1) << '\t'
              << v.reference_bases();
    for (const auto& a : v.alternate_bases()) std::cout << '\t' << a;
    std::cout << '\t' << argmax;
    // alt_allele_indices: which alt-subset this CVO scored.
    std::cout << "\tAAI=";
    for (int i = 0; i < cvo.alt_allele_indices().indices_size(); ++i) {
      if (i) std::cout << ',';
      std::cout << cvo.alt_allele_indices().indices(i);
    }
    // AD and DP from the first call's info, so we can diff variant_caller
    // output between us and upstream at the candidate-emission layer.
    if (v.calls_size() > 0) {
      const auto& info = v.calls(0).info();
      auto it_dp = info.find("DP");
      auto it_ad = info.find("AD");
      std::cout << "\tDP=";
      if (it_dp != info.end() && it_dp->second.values_size() > 0) {
        std::cout << it_dp->second.values(0).int_value();
      }
      std::cout << "\tAD=";
      if (it_ad != info.end()) {
        for (int i = 0; i < it_ad->second.values_size(); ++i) {
          if (i) std::cout << ',';
          std::cout << it_ad->second.values(i).int_value();
        }
      }
    }
    // Append all probabilities at full precision so we can diff against
    // upstream's intermediate CVOs at the postprocess input layer.
    for (int i = 0; i < cvo.genotype_probabilities_size(); ++i) {
      std::cout << '\t' << std::scientific << std::setprecision(17)
                << cvo.genotype_probabilities(i);
    }
    std::cout << '\n';
    ++total;
  }
  std::fprintf(stderr, "%ld records\n", total);
  return 0;
}
