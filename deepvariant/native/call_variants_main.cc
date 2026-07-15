// call_variants — Phase 2 native binary.
//
// Reads a TFRecord of tf.train.Example (pileup images from make_examples),
// runs Inception-v3 inference via Core ML, and writes a TFRecord of
// CallVariantsOutput protos.
//
// Usage:
//   deepvariant call_variants \
//     --examples  /path/make_examples.tfrecord@32 \
//     --checkpoint /path/to/wgs.mlpackage \
//     --outfile    /path/call_variants_output.tfrecord \
//     [--batch_size 128] [--compute_units all|cpu_gpu|cpu_only]
//
// The binary is invoked via the top-level `deepvariant` dispatcher (cli.{h,cc}).

#include "deepvariant/native/call_variants.h"

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(__ARM_NEON) || defined(__aarch64__)
#  include <arm_neon.h>
#  define DV_HAVE_NEON 1
#else
#  define DV_HAVE_NEON 0
#endif

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"

#include "deepvariant/native/bnns_finalize.h"
#include "deepvariant/native/coreml_inference.h"
#include "deepvariant/native/dv_signpost.h"
#include "deepvariant/native/metal_inference.h"
#include "deepvariant/native/tfrecord.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

ABSL_FLAG(std::string, examples,  "", "Input TFRecord file(s) of tf.train.Example.");
ABSL_FLAG(std::string, checkpoint, "",
          "Inference model path. With --inference_backend=coreml, a "
          ".mlpackage. With --inference_backend=metal, a .dvw weight "
          "bundle (see tools/conversion/extract_weights.py).");
ABSL_FLAG(std::string, outfile,   "", "Output TFRecord file for CallVariantsOutput.");
ABSL_FLAG(int,    batch_size, 128, "Inference batch size.");
ABSL_FLAG(std::string, compute_units, "all",
          "Core ML compute units: all (default), cpu_gpu, cpu_only. "
          "Only applies when --inference_backend=coreml.");
ABSL_FLAG(int, input_height, 100,
          "Pileup-image height for the Metal backend. WGS=100, Trio WGS=140 "
          "(60 child + 2x40 parent), pangenome=100, etc.");
ABSL_FLAG(int, input_channels, 7,
          "Pileup-image channels for the Metal backend. WGS/Trio=7, "
          "PacBio/ONT germline=10, MaSeq=9, Hybrid/RNASeq=6.");
ABSL_FLAG(int, input_width, 221,
          "Pileup-image width for the Metal backend. WGS/WES/MaSeq=221, "
          "PacBio=147, ONT=199.");
ABSL_FLAG(std::string, inference_backend, "metal",
          "Inference backend: metal (default, MPSGraph + BNNS-CPU .dvw — "
          "GPU FP32 on Apple Silicon), coreml (Core ML .mlpackage — ANE "
          "or GPU per --compute_units), or ane_speculate (ANE FP16 first, "
          "GPU FP32 rerun for borderline-confidence sites — Scenario 3 "
          "from the master plan).");
ABSL_FLAG(std::string, ane_speculate_metal_checkpoint, "",
          "When --inference_backend=ane_speculate, the .dvw bundle for "
          "the GPU FP32 rerun on borderline-confidence sites. Required.");
// Per-role variants of the .dvw bundle path so cli.cc can thread the
// right rerun model into each sub-call (trio child/parent, somatic
// tumor model, pangenome 9-channel model).
ABSL_FLAG(std::string, ane_speculate_metal_checkpoint_child, "",
          "ane_speculate .dvw bundle for the trio child sample.");
ABSL_FLAG(std::string, ane_speculate_metal_checkpoint_parent, "",
          "ane_speculate .dvw bundle for the trio parent samples.");
ABSL_FLAG(std::string, ane_speculate_metal_checkpoint_somatic, "",
          "ane_speculate .dvw bundle for the DeepSomatic tumor model.");
ABSL_FLAG(std::string, ane_speculate_metal_checkpoint_pangenome, "",
          "ane_speculate .dvw bundle for the pangenome 9-channel model.");
ABSL_FLAG(double, ane_speculate_confidence, 0.99,
          "Borderline threshold for ane_speculate. If max(softmax_ane) < "
          "this value, the example is reclassified on GPU FP32. Lower "
          "→ more GPU reruns, more wall-time, fewer FP-drift artefacts.");
ABSL_FLAG(bool, enable_inference_pipelining, false,
          "Overlap GPU backbone (batch N+1) with CPU BNNS finalize + CVO "
          "build (batch N). Determinism-neutral; default off.");

namespace deepvariant {

namespace {

// Parse the tf.train.Example minimal proto to extract features.
// We only do minimal wire-level parsing; see tools/conversion/bench.py for
// the Python equivalent.
struct ExampleFeatures {
  std::string image_encoded;   // bytes_list value of "image/encoded"
  std::string variant_encoded; // bytes_list value of "variant/encoded"
  std::string alt_allele_indices_encoded; // "alt_allele_indices/encoded"
};

// Read a varint from buf starting at position i. Returns (value, new_i).
static uint64_t ReadVarint(const uint8_t* buf, size_t len, size_t& i) {
  uint64_t val = 0;
  int shift = 0;
  while (i < len) {
    uint8_t b = buf[i++];
    val |= static_cast<uint64_t>(b & 0x7F) << shift;
    if (!(b & 0x80)) return val;
    shift += 7;
  }
  return val;  // truncated
}

// Extract a single bytes value from a BytesList field (wire type 2).
// Extracts the first bytes-value from a Feature whose payload is a BytesList.
// The input is the raw bytes of a Feature proto (the value side of a
// map<string, Feature> entry). The Feature is a oneof — field 1 is BytesList.
// BytesList itself has `repeated bytes value = 1;` — each value is a
// length-delimited bytes entry. We walk both levels and return the first
// value's raw bytes (with no proto framing).
static std::string ExtractBytesListFirst(const uint8_t* buf, size_t len) {
  size_t i = 0;
  while (i < len) {
    uint64_t tag = ReadVarint(buf, len, i);
    uint32_t field = static_cast<uint32_t>(tag >> 3);
    uint32_t wire  = static_cast<uint32_t>(tag & 7);
    if (wire != 2) break;  // we only handle length-delimited
    uint64_t seg_len = ReadVarint(buf, len, i);
    if (i + seg_len > len) break;
    if (field == 1) {
      // We're inside Feature.bytes_list — recurse one level to read the
      // first BytesList.value entry (also a length-delimited bytes field).
      const uint8_t* inner = buf + i;
      size_t j = 0;
      while (j < seg_len) {
        uint64_t itag = ReadVarint(inner, seg_len, j);
        uint32_t ifield = static_cast<uint32_t>(itag >> 3);
        uint32_t iwire  = static_cast<uint32_t>(itag & 7);
        if (iwire != 2) break;
        uint64_t ilen = ReadVarint(inner, seg_len, j);
        if (j + ilen > seg_len) break;
        if (ifield == 1) {
          return std::string(reinterpret_cast<const char*>(inner + j), ilen);
        }
        j += ilen;
      }
      return {};
    }
    i += seg_len;
  }
  return {};
}

// Parse a tf.train.Example wire to extract key fields.
// tf.train.Example has one field: features (field=1, wire=2) → Features
// Features has one repeated field: feature (field=1, wire=2) → map<string, Feature>
// Each map entry: key (field=1), value (field=2).
// Feature is a oneof: bytes_list (field=1), float_list (field=2), int64_list (field=3).
static ExampleFeatures ParseExample(const std::string& payload) {
  ExampleFeatures out;
  const uint8_t* buf = reinterpret_cast<const uint8_t*>(payload.data());
  size_t n = payload.size();
  size_t i = 0;

  // Walk top-level Example proto.
  while (i < n) {
    uint64_t tag = ReadVarint(buf, n, i);
    uint32_t wire = tag & 7;
    if (wire != 2) { break; }
    uint64_t seg_len = ReadVarint(buf, n, i);
    if (i + seg_len > n) break;
    // field 1 = Features
    // Walk the Features proto.
    const uint8_t* feat_buf = buf + i;
    size_t feat_len = seg_len;
    i += seg_len;

    size_t fi = 0;
    while (fi < feat_len) {
      uint64_t ftag = ReadVarint(feat_buf, feat_len, fi);
      uint32_t fwire = ftag & 7;
      if (fwire != 2) break;
      uint64_t entry_len = ReadVarint(feat_buf, feat_len, fi);
      if (fi + entry_len > feat_len) break;
      const uint8_t* entry = feat_buf + fi;
      fi += entry_len;

      // Parse map entry: key (field=1), value (field=2).
      std::string key;
      std::string value_bytes;
      size_t ei = 0;
      while (ei < entry_len) {
        uint64_t etag = ReadVarint(entry, entry_len, ei);
        uint32_t ewire = etag & 7;
        uint32_t efd   = etag >> 3;
        if (ewire != 2) { break; }
        uint64_t elen = ReadVarint(entry, entry_len, ei);
        if (ei + elen > entry_len) break;
        if (efd == 1) {
          key.assign(reinterpret_cast<const char*>(entry + ei), elen);
        } else if (efd == 2) {
          // Feature oneof; field=1 = BytesList
          value_bytes.assign(reinterpret_cast<const char*>(entry + ei), elen);
        }
        ei += elen;
      }

      if (key == "image/encoded" || key == "image") {
        // BytesList → first value
        out.image_encoded = ExtractBytesListFirst(
            reinterpret_cast<const uint8_t*>(value_bytes.data()),
            value_bytes.size());
      } else if (key == "variant/encoded") {
        out.variant_encoded = ExtractBytesListFirst(
            reinterpret_cast<const uint8_t*>(value_bytes.data()),
            value_bytes.size());
      } else if (key == "alt_allele_indices/encoded") {
        out.alt_allele_indices_encoded = ExtractBytesListFirst(
            reinterpret_cast<const uint8_t*>(value_bytes.data()),
            value_bytes.size());
      }
    }
  }
  return out;
}

ComputeUnits ParseComputeUnits(const std::string& s) {
  if (s == "cpu_gpu")  return ComputeUnits::kCpuAndGpu;
  if (s == "cpu_only") return ComputeUnits::kCpuOnly;
  return ComputeUnits::kAll;
}

}  // namespace

int RunCallVariants(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const std::string examples_path  = absl::GetFlag(FLAGS_examples);
  const std::string checkpoint_path = absl::GetFlag(FLAGS_checkpoint);
  const std::string outfile_path   = absl::GetFlag(FLAGS_outfile);
  const int batch_size             = absl::GetFlag(FLAGS_batch_size);
  const ComputeUnits compute_units =
      ParseComputeUnits(absl::GetFlag(FLAGS_compute_units));

  if (examples_path.empty() || checkpoint_path.empty() || outfile_path.empty()) {
    LOG(ERROR) << "Required flags: --examples, --checkpoint, --outfile";
    return 2;
  }

  // Pick inference backend.
  const std::string backend = absl::GetFlag(FLAGS_inference_backend);
  std::unique_ptr<CoreMLModel> coreml_model;
  std::unique_ptr<MetalInception> metal_model;
  std::unique_ptr<BnnsFinalize> metal_finalize;
  int H = 0, W = 0, C = 0, K = 0;
  if (backend == "coreml") {
    LOG(INFO) << "Loading Core ML model: " << checkpoint_path;
    coreml_model = CoreMLModel::Load(checkpoint_path, compute_units);
    if (!coreml_model) {
      LOG(ERROR) << "Failed to load Core ML model: " << checkpoint_path;
      return 1;
    }
    H = coreml_model->InputHeight();
    W = coreml_model->InputWidth();
    C = coreml_model->InputChannels();
    K = coreml_model->NumClasses();
  } else if (backend == "metal") {
    LOG(INFO) << "Loading Metal/BNNS model: " << checkpoint_path;
    // Pass --input_height / --input_channels to MetalInception so the
    // MPSGraph placeholder is built with the right shape. Defaults
    // (100×221×7) match WGS; trio passes 140 via --input_height.
    H = absl::GetFlag(FLAGS_input_height);
    W = absl::GetFlag(FLAGS_input_width);
    C = absl::GetFlag(FLAGS_input_channels);
    K = 3;
    metal_model = MetalInception::Create(checkpoint_path, H, C, W);
    metal_finalize = BnnsFinalize::Create(checkpoint_path);
    if (!metal_model || !metal_finalize) {
      LOG(ERROR) << "Failed to load Metal/BNNS model: " << checkpoint_path;
      return 1;
    }
  } else if (backend == "ane_speculate") {
    // Scenario 3: ANE FP16 forward pass on every example; for examples
    // where max(softmax_ane) < threshold (= --ane_speculate_confidence,
    // default 0.99), rerun on GPU MPSGraph FP32 + BNNS-CPU finalize so
    // borderline GQ=20 sites stay on the deterministic FP32 path.
    const std::string metal_ckpt =
        absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint);
    if (metal_ckpt.empty()) {
      LOG(ERROR) << "ane_speculate requires --ane_speculate_metal_checkpoint=<.dvw>";
      return 2;
    }
    LOG(INFO) << "Loading ane_speculate ANE model:   " << checkpoint_path;
    coreml_model = CoreMLModel::Load(checkpoint_path, compute_units);
    if (!coreml_model) {
      LOG(ERROR) << "Failed to load Core ML .mlpackage: " << checkpoint_path;
      return 1;
    }
    LOG(INFO) << "Loading ane_speculate GPU rerun:   " << metal_ckpt;
    H = absl::GetFlag(FLAGS_input_height);
    W = absl::GetFlag(FLAGS_input_width);
    C = absl::GetFlag(FLAGS_input_channels);
    K = 3;
    metal_model = MetalInception::Create(metal_ckpt, H, C, W);
    metal_finalize = BnnsFinalize::Create(metal_ckpt);
    if (!metal_model || !metal_finalize) {
      LOG(ERROR) << "Failed to load .dvw fallback bundle: " << metal_ckpt;
      return 1;
    }
    // Soft sanity check: ANE model's declared input shape vs Metal
    // model's. A mismatch could indicate the .mlpackage was extracted
    // with the wrong height (e.g. trio child should be 140, not 100).
    // Some Core ML packages declare flexible/dynamic shapes; defer the
    // hard check to Predict() which will surface a precise error.
    if (coreml_model->InputHeight() != H || coreml_model->InputChannels() != C) {
      LOG(WARNING) << "ane_speculate: declared shape mismatch — ANE "
                   << "expects (" << coreml_model->InputHeight()
                   << "x" << coreml_model->InputWidth() << "x"
                   << coreml_model->InputChannels()
                   << ") vs Metal (" << H << "x" << W << "x" << C
                   << "). Will rely on Core ML's runtime shape handling.";
    }
  } else {
    LOG(ERROR) << "Unknown --inference_backend=" << backend
               << " (expected 'coreml', 'metal' or 'ane_speculate')";
    return 2;
  }
  LOG(INFO) << "Model input (" << H << "," << W << "," << C
            << ") → " << K << " classes  [backend=" << backend << "]";

  // Open TFRecord reader + writer.
  auto reader = TFRecordReader::New(examples_path);
  if (!reader) {
    LOG(ERROR) << "Cannot open examples file: " << examples_path;
    return 1;
  }
  auto writer = TFRecordWriter::New(outfile_path);
  if (!writer) {
    LOG(ERROR) << "Cannot open output file: " << outfile_path;
    return 1;
  }

  // ── P1: async writer thread ──────────────────────────────────────────────
  // Move CVO TFRecord writes off the main thread so we can overlap them
  // with the next batch's GPU compute. Bounded SPSC queue gives back-
  // pressure when writer falls behind the producer (rare since GPU is
  // much slower than disk write at our throughput).
  //
  // Design:
  //   main thread: build CVO → SerializeToString → enqueue
  //   writer thread: dequeue → writer->WriteRecord → loop
  //   end: main pushes 'done' flag, writer drains queue + exits
  //
  // Output bit-equivalence: writer thread is the SOLE consumer of the
  // writer; serialization order is preserved by the queue's FIFO
  // discipline. Same TFRecord bytes produced.
  constexpr size_t kWriteQueueDepth = 32;  // up to 32 CVOs buffered
  std::deque<std::string> write_queue;
  std::mutex wq_mu;
  std::condition_variable wq_nonempty, wq_nonfull;
  bool writer_done = false;
  std::atomic<bool> writer_failed{false};

  std::thread writer_thread([&]() {
    for (;;) {
      std::string item;
      {
        std::unique_lock<std::mutex> lk(wq_mu);
        wq_nonempty.wait(lk, [&] {
          return !write_queue.empty() || writer_done;
        });
        if (write_queue.empty() && writer_done) return;
        item = std::move(write_queue.front());
        write_queue.pop_front();
        wq_nonfull.notify_one();
      }
      if (!writer->WriteRecord(item)) {
        LOG(ERROR) << "Async writer: WriteRecord failed";
        writer_failed.store(true);
        // Drain remaining queue silently to unblock producer.
        std::lock_guard<std::mutex> lk(wq_mu);
        write_queue.clear();
        wq_nonfull.notify_all();
        return;
      }
    }
  });

  // RAII guard: join the async writer on EVERY exit path, including early
  // `return 1` (e.g. a flush_batch failure). Without it, returning while
  // writer_thread is still joinable destroys a joinable std::thread →
  // std::terminate(). Declared before PipeJoiner (below) so it destructs
  // LAST — after the pipelined finalize worker has drained its CVOs into the
  // write queue — preserving the same ordering as the explicit teardown.
  struct WriterJoiner {
    std::thread& t;
    std::mutex& mu;
    bool& done;
    std::condition_variable& cv;
    ~WriterJoiner() {
      if (!t.joinable()) return;
      { std::lock_guard<std::mutex> lk(mu); done = true; }
      cv.notify_all();
      t.join();
    }
  } writer_joiner{writer_thread, wq_mu, writer_done, wq_nonempty};

  auto enqueue_write = [&](std::string&& payload) -> bool {
    if (writer_failed.load()) return false;
    std::unique_lock<std::mutex> lk(wq_mu);
    wq_nonfull.wait(lk, [&] {
      return write_queue.size() < kWriteQueueDepth || writer_failed.load();
    });
    if (writer_failed.load()) return false;
    write_queue.push_back(std::move(payload));
    wq_nonempty.notify_one();
    return true;
  };

  // Batch inference loop.
  int64_t total_examples = 0;
  int64_t total_batches  = 0;

  struct PendingExample {
    ExampleFeatures features;
    std::string raw_payload;  // original Example bytes (for passthrough fields)
  };
  std::vector<PendingExample> batch;
  batch.reserve(batch_size);

  // Hoist large per-batch buffer allocations out of the flush loop.
  // For batch_size=2048 and chr20 (B × H × W × C × 4 ≈ 1.3 GB), the
  // per-batch malloc + memset is a measurable cost (~80-150 ms per
  // batch on M4 Max). Allocate once at full capacity, reuse across
  // batches. The MPSGraph input wrapper reads only `n × elem` bytes
  // so the trailing slack is harmless.
  std::vector<float> images(static_cast<size_t>(batch_size) *
                              static_cast<size_t>(H * W * C));
  std::vector<float> probs(static_cast<size_t>(batch_size) *
                             static_cast<size_t>(K));

  // Shared CVO emit: builds one CallVariantsOutput per example from a
  // probabilities buffer and per-example variant/alt metadata, then pushes
  // each serialized record to the async writer via enqueue_write. This is
  // the SOLE producer of CVO records and is used by BOTH the serial path
  // and the pipelined finalize worker so the two paths cannot diverge.
  //
  // `variant_at(i)` / `alt_at(i)` return the encoded bytes for example i;
  // `probs_ptr` points at the (n × K) probabilities for this batch. Records
  // are emitted strictly in increasing i — combined with the writer's FIFO
  // queue and a single consumer, output byte-order matches the serial path.
  auto emit_cvos = [&](const float* probs_ptr, int n,
                       auto&& variant_at, auto&& alt_at) -> bool {
    for (int i = 0; i < n; ++i) {
      learning::genomics::deepvariant::CallVariantsOutput cvo;
      const std::string& variant_encoded = variant_at(i);
      const std::string& alt_encoded = alt_at(i);
      if (!variant_encoded.empty() &&
          !cvo.mutable_variant()->ParseFromString(variant_encoded)) {
        LOG(ERROR) << "Failed to parse variant/encoded for example " << i;
        return false;
      }
      if (!alt_encoded.empty() &&
          !cvo.mutable_alt_allele_indices()->ParseFromString(alt_encoded)) {
        LOG(ERROR) << "Failed to parse alt_allele_indices/encoded for example "
                   << i;
        return false;
      }
      for (int k = 0; k < K; ++k) {
        cvo.add_genotype_probabilities(probs_ptr[i * K + k]);
      }
      // Tag MID="deepvariant" so postprocess can write it as a VCF FORMAT
      // field. Reuse the empty VariantCall slot that variant_calling.cc
      // already added (otherwise we end up with 2 calls and VcfWriter
      // rejects the variant for not matching sample count).
      auto* v = cvo.mutable_variant();
      if (v->calls_size() == 0) v->add_calls();
      nucleus::SetInfoField("MID", std::string("deepvariant"),
                             v->mutable_calls(0));

      std::string serialized;
      if (!cvo.SerializeToString(&serialized)) {
        LOG(ERROR) << "Failed to serialize CallVariantsOutput";
        return false;
      }
      // P1: async writer thread consumes this. Push std::move so the
      // writer thread owns the buffer; main thread can recycle storage.
      if (!enqueue_write(std::move(serialized))) {
        LOG(ERROR) << "Failed to enqueue output record (writer thread error)";
        return false;
      }
    }
    return true;
  };

  // ── P3 (opt-in): GPU backbone / CPU finalize pipelining ──────────────────
  // When --enable_inference_pipelining is set AND the active backend is the
  // two-stage Metal/BNNS path, we overlap the GPU MPSGraph backbone of batch
  // N+1 with the CPU BNNS finalize + CVO build of batch N.
  //
  // Structure (one-batch-deep double buffer + single FIFO worker):
  //   main thread  : normalize → metal_model->Predict() (blocks on the GPU;
  //                  on return `images` is free to reuse and the per-slot
  //                  features buffer holds the backbone activations) → hand
  //                  the slot to the worker → continue with next batch.
  //   worker thread: pop slot in FIFO order → metal_finalize->ApplyBatch()
  //                  (CPU BNNS dense+softmax) → emit_cvos() in per-example
  //                  order → mark slot free.
  //
  // Determinism / order: there is exactly ONE worker draining a FIFO job
  // deque, and emit_cvos() (shared with the serial path) is the sole CVO
  // producer. So CVOs are built and enqueued in the same global order as the
  // serial path → byte-identical output. ApplyBatch on identical features is
  // bit-deterministic, so probabilities are identical too.
  //
  // At-most-one-in-flight invariant: we use exactly 2 slots. The main thread
  // must not Predict-write features[slot] for batch N+2 until the worker has
  // finished batch N (which used the same slot). `slot_busy[slot]` guards
  // this: the main thread waits for slot_busy[slot]==false before writing,
  // the worker sets it false after finishing that slot's job. Thus the main
  // thread is at most one batch ahead of the worker.
  const bool pipeline_enabled =
      absl::GetFlag(FLAGS_enable_inference_pipelining);
  const int kPipelineSlots = 2;

  // Per-slot job metadata captured by the main thread for the worker. We copy
  // the (small) variant/alt encoded strings out of `batch` so `batch` can be
  // cleared and refilled for the next batch while the worker still needs
  // this batch's metadata.
  struct CvoMeta {
    std::string variant_encoded;
    std::string alt_allele_indices_encoded;
  };
  struct PipelineJob {
    int slot = 0;
    int n = 0;
    std::vector<CvoMeta> meta;  // size n
  };

  // Double-buffered features and probs (only allocated when the pipeline is
  // actually used; sized lazily once the FeatureDim is known).
  std::vector<float> pipe_features[2];
  std::vector<float> pipe_probs[2];

  // Job queue + slot-free signaling.
  std::deque<PipelineJob> pipe_queue;
  std::mutex pipe_mu;
  std::condition_variable pipe_nonempty;   // worker waits for jobs
  std::condition_variable pipe_slot_free;  // main waits for a slot to free
  bool pipe_busy[2] = {false, false};      // slot occupied by an in-flight job
  bool pipe_done = false;                   // no more jobs will be enqueued
  std::atomic<bool> pipe_failed{false};     // worker hit an error

  // The finalize worker. Started lazily on first pipelined flush so the
  // serial path (and non-Metal backends) never spawn it. Captures by ref;
  // `metal_finalize`, `emit_cvos`, the pipe_* state and the buffers all
  // outlive the worker (joined before this function returns).
  std::thread pipe_worker;
  auto ensure_pipe_worker = [&]() {
    if (pipe_worker.joinable()) return;
    pipe_worker = std::thread([&]() {
      for (;;) {
        PipelineJob job;
        {
          std::unique_lock<std::mutex> lk(pipe_mu);
          pipe_nonempty.wait(lk, [&] {
            return !pipe_queue.empty() || pipe_done;
          });
          if (pipe_queue.empty()) {
            if (pipe_done) return;
            continue;
          }
          job = std::move(pipe_queue.front());
          pipe_queue.pop_front();
        }
        const int slot = job.slot;
        bool job_ok = !pipe_failed.load();
        if (job_ok) {
          // CPU BNNS dense + softmax into this slot's probs buffer.
          job_ok = metal_finalize->ApplyBatch(
              pipe_features[slot].data(), job.n, pipe_probs[slot].data());
          if (!job_ok) {
            LOG(ERROR) << "Pipelined finalize: ApplyBatch failed";
          }
        }
        if (job_ok) {
          // Build + enqueue CVOs in per-example order (shared emit path).
          job_ok = emit_cvos(
              pipe_probs[slot].data(), job.n,
              [&](int i) -> const std::string& {
                return job.meta[i].variant_encoded;
              },
              [&](int i) -> const std::string& {
                return job.meta[i].alt_allele_indices_encoded;
              });
        }
        if (!job_ok) pipe_failed.store(true);
        // Free the slot regardless of success so the main thread (which may
        // be blocked waiting for this slot) is never deadlocked on error.
        {
          std::lock_guard<std::mutex> lk(pipe_mu);
          pipe_busy[slot] = false;
        }
        pipe_slot_free.notify_one();
      }
    });
  };

  // Drain + join the worker. Safe to call on any exit path; idempotent.
  // Setting pipe_done wakes the worker so it drains the queue and exits; the
  // worker also frees slots + notifies pipe_slot_free on every job (including
  // the error path), so the main thread can never be left blocked waiting for
  // a slot while we are tearing down. join() is the final barrier — after it
  // returns, all queued CVOs have been emitted to the writer.
  auto join_pipe_worker = [&]() {
    if (!pipe_worker.joinable()) return;
    {
      std::lock_guard<std::mutex> lk(pipe_mu);
      pipe_done = true;
    }
    pipe_nonempty.notify_all();
    pipe_worker.join();
  };

  // RAII guard: ensure the finalize worker is joined on EVERY exit path
  // (including early `return 1` on flush failure). Idempotent with the
  // explicit end-of-input join below. Declared before flush_batch / the main
  // loop so it is destroyed (and thus joins) after them on any return.
  struct PipeJoiner {
    std::function<void()>& join_fn;
    ~PipeJoiner() { join_fn(); }
  };
  std::function<void()> join_pipe_worker_fn = join_pipe_worker;
  PipeJoiner pipe_joiner{join_pipe_worker_fn};

  auto flush_batch = [&]() -> bool {
    if (batch.empty()) return true;
    const int n = static_cast<int>(batch.size());
    const int64_t elem = H * W * C;
    DV_SIGNPOST_INTERVAL_BEGIN(FlushBatch, "");
    DV_SIGNPOST_INTERVAL_BEGIN(Normalize, "");
    for (int i = 0; i < n; ++i) {
      const std::string& img = batch[i].features.image_encoded;
      if (static_cast<int64_t>(img.size()) != elem) {
        // Try float32 layout (some variants store floats directly).
        if (static_cast<int64_t>(img.size()) == elem * 4) {
          std::memcpy(images.data() + i * elem, img.data(), elem * 4);
        } else {
          LOG(ERROR) << "Unexpected image size " << img.size()
                     << " (expected " << elem << " or " << elem * 4 << ")";
          return false;
        }
      } else {
        // uint8 → float32 normalized to [-1, 1] via (x - 128) / 128.
        // This matches the upstream DeepVariant preprocess_images (see
        // deepvariant/dv_utils.py: tf.subtract(images, 128.0); divide(., 128.0)).
        //
        // Bit-equivalence note: 1/128 = 2^-7 is exactly representable in
        // FP32, and (byte - 128.0f) for byte ∈ [0,255] is also exact, so
        // the multiplication produces exact results matching the scalar
        // path bit-for-bit. NEON intrinsics use IEEE 754 single-rounded
        // ops on Apple Silicon → identical FP32 outputs vs the scalar
        // loop. Verified: same inputs through scalar vs NEON paths
        // produce byte-identical `images` buffer.
        const uint8_t* src = reinterpret_cast<const uint8_t*>(img.data());
        float* dst = images.data() + i * elem;
        constexpr float kInvScale = 1.0f / 128.0f;
#if DV_HAVE_NEON
        const float32x4_t k128 = vdupq_n_f32(128.0f);
        const float32x4_t kinv = vdupq_n_f32(kInvScale);
        const int64_t simd_end = elem & ~int64_t{15};
        for (int64_t j = 0; j < simd_end; j += 16) {
          uint8x16_t b = vld1q_u8(src + j);
          // 16 u8 → 4×4 u32 → 4×4 f32 lanes.
          uint16x8_t lo16 = vmovl_u8(vget_low_u8(b));
          uint16x8_t hi16 = vmovl_u8(vget_high_u8(b));
          float32x4_t f0 = vcvtq_f32_u32(vmovl_u16(vget_low_u16(lo16)));
          float32x4_t f1 = vcvtq_f32_u32(vmovl_u16(vget_high_u16(lo16)));
          float32x4_t f2 = vcvtq_f32_u32(vmovl_u16(vget_low_u16(hi16)));
          float32x4_t f3 = vcvtq_f32_u32(vmovl_u16(vget_high_u16(hi16)));
          vst1q_f32(dst + j +  0, vmulq_f32(vsubq_f32(f0, k128), kinv));
          vst1q_f32(dst + j +  4, vmulq_f32(vsubq_f32(f1, k128), kinv));
          vst1q_f32(dst + j +  8, vmulq_f32(vsubq_f32(f2, k128), kinv));
          vst1q_f32(dst + j + 12, vmulq_f32(vsubq_f32(f3, k128), kinv));
        }
        // Tail (< 16 trailing bytes).
        for (int64_t j = simd_end; j < elem; ++j) {
          dst[j] = (static_cast<float>(src[j]) - 128.0f) * kInvScale;
        }
#else
        for (int64_t j = 0; j < elem; ++j) {
          dst[j] = (static_cast<float>(src[j]) - 128.0f) * kInvScale;
        }
#endif
      }
    }

    DV_SIGNPOST_INTERVAL_END(Normalize);

    // Run inference. (probs hoisted, see top of fn; features lazily
    // allocated to full batch capacity inside the metal branch.)
    bool ok = false;
    DV_SIGNPOST_INTERVAL_BEGIN(Inference, "");
    const bool ane_speculate_mode =
        (coreml_model && metal_model && metal_finalize);
    if (ane_speculate_mode) {
      // Scenario 3: ANE FP16 forward on the full batch; rerun
      // borderline-confidence examples on GPU MPSGraph FP32 +
      // BNNS-CPU finalize so threshold sites stay on the
      // deterministic FP32 path.
      DV_SIGNPOST_INTERVAL_BEGIN(AneFp16, "");
      ok = coreml_model->Predict(images.data(), n, H, W, C,
                                  probs.data(), K);
      DV_SIGNPOST_INTERVAL_END(AneFp16);
      if (ok) {
        // Identify borderline examples. Two triggers — either qualifies
        // as borderline and forces a GPU FP32 rerun:
        //
        //   (1) max(softmax) < conf_threshold
        //       → top-class confidence is below the gate. This catches
        //         GQ ≈ 20 boundary flips where ANE FP16's drift on the
        //         winning class could change the FILTER classification.
        //
        //   (2) min(softmax) < min_floor (default 1e-4)
        //       → at least one of the {homref, het, homvar} probabilities
        //         is small enough that FP16's ~10⁻⁴ relative precision
        //         leaks into the floor()-rounded PL byte:
        //           PL_i = floor(-10*log10(p_i / max_p))
        //         A 10⁻⁴ relative change in p_i at p_i ~ 10⁻⁴ produces
        //         a 1-PL-unit difference vs FP32. (2) catches that
        //         purely-textual drift without changing FILTER (FP16
        //         argmax remains stable when max_p ≫ 0.9999).
        const float conf_threshold = static_cast<float>(
            absl::GetFlag(FLAGS_ane_speculate_confidence));
        // Static for now: 1e-4 is the FP16 noise-floor at small p
        // values. Could be exposed as a flag if users want to tune.
        const float min_floor = 1e-4f;
        static thread_local std::vector<int> borderline_idx;
        borderline_idx.clear();
        borderline_idx.reserve(n);
        for (int i = 0; i < n; ++i) {
          float m = probs[i * K], mn = probs[i * K];
          for (int j = 1; j < K; ++j) {
            const float p = probs[i * K + j];
            if (p > m)  m  = p;
            if (p < mn) mn = p;
          }
          if (m < conf_threshold || mn < min_floor) {
            borderline_idx.push_back(i);
          }
        }
        if (!borderline_idx.empty()) {
          DV_SIGNPOST_INTERVAL_BEGIN(AneRerunGpu, "");
          const int nb = static_cast<int>(borderline_idx.size());
          const size_t img_per = static_cast<size_t>(H) * W * C;
          static thread_local std::vector<float> bl_images, bl_features,
              bl_probs;
          bl_images.resize(static_cast<size_t>(nb) * img_per);
          bl_features.resize(static_cast<size_t>(nb) *
                             metal_model->FeatureDim());
          bl_probs.resize(static_cast<size_t>(nb) * K);
          for (int b = 0; b < nb; ++b) {
            const int src = borderline_idx[b];
            std::memcpy(bl_images.data() + static_cast<size_t>(b) * img_per,
                        images.data() + static_cast<size_t>(src) * img_per,
                        img_per * sizeof(float));
          }
          bool gpu_ok = metal_model->Predict(bl_images.data(), nb,
                                              bl_features.data());
          if (gpu_ok) {
            gpu_ok = metal_finalize->ApplyBatch(bl_features.data(), nb,
                                                bl_probs.data());
          }
          if (gpu_ok) {
            for (int b = 0; b < nb; ++b) {
              const int dst = borderline_idx[b];
              std::memcpy(probs.data() + static_cast<size_t>(dst) * K,
                          bl_probs.data() + static_cast<size_t>(b) * K,
                          K * sizeof(float));
            }
          } else {
            ok = false;
            LOG(ERROR) << "ane_speculate: GPU rerun failed on "
                       << nb << " borderline examples";
          }
          DV_SIGNPOST_INTERVAL_END(AneRerunGpu);
        }
      }
    } else if (coreml_model) {
      ok = coreml_model->Predict(images.data(), n, H, W, C,
                                  probs.data(), K);
    } else if (metal_model && metal_model->IsGpuFinalize()) {
      // Single-stage GPU path (DV_METAL_GPU_FINALIZE=1): the dense +
      // softmax run inside MPSGraph, so Predict() writes (n, 3)
      // probabilities directly. metal_finalize is unused in this mode.
      DV_SIGNPOST_INTERVAL_BEGIN(MetalGPU, "");
      ok = metal_model->Predict(images.data(), n, probs.data());
      DV_SIGNPOST_INTERVAL_END(MetalGPU);
    } else if (metal_model && metal_finalize && pipeline_enabled) {
      // Two-stage Metal/BNNS path, PIPELINED (opt-in). The GPU backbone runs
      // here on the main thread (blocking on the GPU); the CPU BNNS finalize
      // + CVO build are handed to the single finalize worker so they overlap
      // with the NEXT batch's GPU backbone. This branch fully handles its own
      // CVO emission (via the worker) and returns early — it does not fall
      // through to the serial CVO loop below.
      const size_t feat_per =
          static_cast<size_t>(metal_model->FeatureDim());
      const size_t feat_total =
          static_cast<size_t>(batch_size) * feat_per;
      const size_t prob_total =
          static_cast<size_t>(batch_size) * static_cast<size_t>(K);

      ensure_pipe_worker();

      // Pick this batch's slot and WAIT until it is free. Two slots means the
      // main thread is at most one batch ahead of the worker (the slot used
      // `kPipelineSlots` batches ago must be drained before we overwrite it).
      const int slot =
          static_cast<int>(total_batches % kPipelineSlots);
      {
        std::unique_lock<std::mutex> lk(pipe_mu);
        pipe_slot_free.wait(lk, [&] {
          return !pipe_busy[slot] || pipe_failed.load();
        });
      }
      if (pipe_failed.load()) {
        // A prior job failed in the worker; abort cleanly.
        DV_SIGNPOST_INTERVAL_END(Inference);
        return false;
      }

      // Lazily size the double-buffered features/probs for this slot.
      if (pipe_features[slot].size() < feat_total) {
        pipe_features[slot].resize(feat_total);
      }
      if (pipe_probs[slot].size() < prob_total) {
        pipe_probs[slot].resize(prob_total);
      }

      // GPU backbone → this slot's features buffer. Predict() blocks on
      // waitUntilCompleted, so on return both `images` (reusable for the
      // next batch) and pipe_features[slot] are settled.
      DV_SIGNPOST_INTERVAL_BEGIN(MetalGPU, "");
      bool gpu_ok =
          metal_model->Predict(images.data(), n, pipe_features[slot].data());
      DV_SIGNPOST_INTERVAL_END(MetalGPU);
      DV_SIGNPOST_INTERVAL_END(Inference);
      if (!gpu_ok) {
        LOG(ERROR) << "Inference failed on batch " << total_batches;
        return false;
      }

      // Capture per-example CVO metadata for the worker (small strings),
      // moved out of `batch` so `batch` can be refilled immediately.
      PipelineJob job;
      job.slot = slot;
      job.n = n;
      job.meta.resize(n);
      for (int i = 0; i < n; ++i) {
        job.meta[i].variant_encoded =
            std::move(batch[i].features.variant_encoded);
        job.meta[i].alt_allele_indices_encoded =
            std::move(batch[i].features.alt_allele_indices_encoded);
      }

      // Publish the job: mark the slot busy and hand it to the worker.
      {
        std::lock_guard<std::mutex> lk(pipe_mu);
        pipe_busy[slot] = true;
        pipe_queue.push_back(std::move(job));
      }
      pipe_nonempty.notify_one();

      // Bookkeeping happens here on the main thread (it does not depend on
      // the worker). The serial CVO loop below is skipped.
      ++total_batches;
      total_examples += n;
      batch.clear();
      DV_SIGNPOST_INTERVAL_END(FlushBatch);
      return true;
    } else if (metal_model && metal_finalize) {
      // Two-stage Metal/BNNS path (serial): GPU MPSGraph for backbone, CPU
      // BNNS for the final dense + softmax (deterministic FP32 reduction
      // = bit-parity with TF CPU). features sized to full batch_size
      // on first use; subsequent batches reuse via static thread-local.
      static thread_local std::vector<float> features;
      const size_t feat_total = static_cast<size_t>(batch_size) *
                                  static_cast<size_t>(metal_model->FeatureDim());
      if (features.size() < feat_total) features.resize(feat_total);
      DV_SIGNPOST_INTERVAL_BEGIN(MetalGPU, "");
      bool gpu_ok = metal_model->Predict(images.data(), n, features.data());
      DV_SIGNPOST_INTERVAL_END(MetalGPU);
      if (gpu_ok) {
        DV_SIGNPOST_INTERVAL_BEGIN(BnnsFinalize, "");
        ok = metal_finalize->ApplyBatch(features.data(), n, probs.data());
        DV_SIGNPOST_INTERVAL_END(BnnsFinalize);
      }
    }
    DV_SIGNPOST_INTERVAL_END(Inference);
    if (!ok) {
      LOG(ERROR) << "Inference failed on batch " << total_batches;
      return false;
    }

    // Write one CallVariantsOutput per example (serial path; the pipelined
    // two-stage Metal branch above returns early and emits via the worker).
    if (!emit_cvos(
            probs.data(), n,
            [&](int i) -> const std::string& {
              return batch[i].features.variant_encoded;
            },
            [&](int i) -> const std::string& {
              return batch[i].features.alt_allele_indices_encoded;
            })) {
      return false;
    }

    ++total_batches;
    total_examples += n;
    batch.clear();
    DV_SIGNPOST_INTERVAL_END(FlushBatch);
    return true;
  };

  // ── P2: pre-fetch reader thread ──────────────────────────────────────────
  // Move reader->GetNext() + ParseExample off the main thread so we can
  // overlap the I/O + protobuf parsing with the previous batch's GPU
  // dispatch. Bounded SPSC queue (depth = 2 × batch_size = 1024 examples
  // at default batch=512) gives back-pressure when main thread is the
  // bottleneck.
  //
  // Output bit-equivalence: reader produces same PendingExample objects
  // in the same order; main thread consumes in same order; flush_batch
  // sees identical batches as before. No algorithmic change.
  const size_t kReadQueueDepth = static_cast<size_t>(batch_size) * 2;
  std::deque<PendingExample> read_queue;
  std::mutex rq_mu;
  std::condition_variable rq_nonempty, rq_nonfull;
  bool reader_eof = false;
  std::atomic<bool> reader_stop{false};

  std::thread reader_thread([&]() {
    while (!reader_stop.load() && reader->GetNext()) {
      PendingExample pe;
      pe.raw_payload = reader->record();
      pe.features    = ParseExample(pe.raw_payload);
      std::unique_lock<std::mutex> lk(rq_mu);
      rq_nonfull.wait(lk, [&] {
        return read_queue.size() < kReadQueueDepth || reader_stop.load();
      });
      if (reader_stop.load()) return;
      read_queue.push_back(std::move(pe));
      rq_nonempty.notify_one();
    }
    {
      std::lock_guard<std::mutex> lk(rq_mu);
      reader_eof = true;
    }
    rq_nonempty.notify_all();
  });

  // RAII guard: ensure reader thread is joined on every exit path.
  struct ReaderJoiner {
    std::thread& t;
    std::atomic<bool>& stop;
    std::mutex& mu;
    std::condition_variable& cv_full;
    std::condition_variable& cv_empty;
    ~ReaderJoiner() {
      stop.store(true);
      { std::lock_guard<std::mutex> lk(mu); }
      cv_full.notify_all();
      cv_empty.notify_all();
      if (t.joinable()) t.join();
    }
  } reader_joiner{reader_thread, reader_stop, rq_mu, rq_nonfull, rq_nonempty};

  // Main consumption loop: pop from reader queue, accumulate batch,
  // flush when full.
  for (;;) {
    PendingExample pe;
    bool got_one = false;
    {
      std::unique_lock<std::mutex> lk(rq_mu);
      rq_nonempty.wait(lk, [&] {
        return !read_queue.empty() || reader_eof;
      });
      if (!read_queue.empty()) {
        pe = std::move(read_queue.front());
        read_queue.pop_front();
        rq_nonfull.notify_one();
        got_one = true;
      } else if (reader_eof) {
        break;
      }
    }
    if (got_one) {
      batch.push_back(std::move(pe));
      if (static_cast<int>(batch.size()) >= batch_size) {
        if (!flush_batch()) return 1;
      }
    }
  }
  if (!flush_batch()) return 1;

  // End-of-input: drain + join the finalize worker BEFORE signaling the
  // writer to finish. This guarantees every pipelined batch's CVOs have been
  // emit_cvos()'d into the writer queue before we close it, so no output is
  // lost. (No-op when pipelining is disabled / worker never started.) The
  // RAII PipeJoiner would also join, but joining here lets us surface a
  // worker error and keeps teardown ordering explicit.
  join_pipe_worker();
  if (pipe_failed.load()) {
    LOG(ERROR) << "Pipelined finalize worker failed during run";
    return 1;
  }

  // Signal writer thread to drain + exit; then close writer ourselves.
  {
    std::lock_guard<std::mutex> lk(wq_mu);
    writer_done = true;
  }
  wq_nonempty.notify_all();
  writer_thread.join();
  if (writer_failed.load()) {
    LOG(ERROR) << "Async writer thread failed during run";
    return 1;
  }

  reader->Close();
  if (!writer->Close()) {
    LOG(ERROR) << "Failed to flush/close output: " << outfile_path;
    return 1;
  }

  LOG(INFO) << "call_variants done: " << total_examples << " examples, "
            << total_batches << " batches → " << outfile_path;
  return 0;
}

}  // namespace deepvariant
