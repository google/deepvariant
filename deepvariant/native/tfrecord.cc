// TFRecord reader/writer for deepvariant native runtime.
// See tfrecord.h for the format description.

#include "deepvariant/native/tfrecord.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>

#include "absl/crc/crc32c.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"

namespace deepvariant {

namespace {
constexpr uint32_t kMaskDelta = 0xa282ead8UL;

uint32_t MaskedCrc32c(const char* data, size_t n) {
  uint32_t crc = static_cast<uint32_t>(absl::ComputeCrc32c({data, n}));
  return ((crc >> 15) | (crc << 17)) + kMaskDelta;
}

// Expand `prefix@N` to {prefix-00000-of-NNNNN, ..., prefix-(N-1)-of-NNNNN}.
// Plain paths (no `@`) pass through as a single-element list.
std::vector<std::string> ExpandShards(const std::string& spec) {
  auto at = spec.find('@');
  if (at == std::string::npos) return {spec};
  const std::string prefix = spec.substr(0, at);
  int n = 0;
  if (!absl::SimpleAtoi(spec.substr(at + 1), &n) || n <= 0) return {spec};
  std::vector<std::string> paths;
  paths.reserve(n);
  for (int i = 0; i < n; ++i) {
    paths.push_back(absl::StrCat(prefix, "-",
                                  absl::Dec(i, absl::kZeroPad5),
                                  "-of-", absl::Dec(n, absl::kZeroPad5)));
  }
  return paths;
}
}  // namespace

std::string ShardName(const std::string& spec, int task_id) {
  auto at = spec.find('@');
  if (at == std::string::npos) return spec;
  const std::string prefix = spec.substr(0, at);
  int n = 0;
  if (!absl::SimpleAtoi(spec.substr(at + 1), &n) || n <= 0) return spec;
  return absl::StrCat(prefix, "-", absl::Dec(task_id, absl::kZeroPad5),
                       "-of-", absl::Dec(n, absl::kZeroPad5));
}

// ---------------------------------------------------------------------------
// TFRecordReader
// ---------------------------------------------------------------------------

// Large streambuf for the reader. The writer coalesces into a 1 MiB buffer
// (kBufBytes), but the reader otherwise uses the default ~4-8 KiB streambuf
// and issues several stream ops per record. A 1 MiB read buffer cuts the
// number of underlying read() syscalls dramatically. Pure buffering: the
// on-disk format and decoded bytes are unchanged.
namespace {
constexpr size_t kReaderBufBytes = 1 << 20;  // 1 MiB read buffer
}  // namespace

struct TFRecordReader::Impl {
  std::vector<std::string> paths;
  size_t current_index = 0;
  // Persistent buffer backing the stream's streambuf. Its lifetime must
  // cover all reads, so it is owned by the Impl alongside the stream.
  std::vector<char> read_buf;
  std::ifstream stream;

  explicit Impl(const std::string& spec)
      : paths(ExpandShards(spec)), read_buf(kReaderBufBytes) {
    if (!paths.empty()) {
      // pubsetbuf must be called BEFORE open to take effect.
      stream.rdbuf()->pubsetbuf(read_buf.data(), read_buf.size());
      stream.open(paths[0], std::ios::binary);
    }
  }

  // Advance to the next shard if the current one is exhausted; returns true
  // if a stream is currently open and ready for reading.
  bool EnsureOpen() {
    if (stream.is_open() && stream.good()) return true;
    while (current_index + 1 < paths.size()) {
      stream.close();
      ++current_index;
      stream.clear();
      // pubsetbuf must be called BEFORE open to take effect.
      stream.rdbuf()->pubsetbuf(read_buf.data(), read_buf.size());
      stream.open(paths[current_index], std::ios::binary);
      if (stream.is_open() && stream.good()) return true;
    }
    return false;
  }
};

TFRecordReader::TFRecordReader() = default;
TFRecordReader::~TFRecordReader() = default;

std::unique_ptr<TFRecordReader> TFRecordReader::New(
    const std::string& path, const std::string& /*compression_type*/) {
  auto impl = std::make_unique<Impl>(path);
  if (impl->paths.empty()) return nullptr;
  if (!impl->stream.is_open()) return nullptr;
  auto r = std::unique_ptr<TFRecordReader>(new TFRecordReader());
  r->impl_ = std::move(impl);
  return r;
}

bool TFRecordReader::GetNext() {
  if (!impl_) return false;
  while (true) {
    auto& s = impl_->stream;
    if (s.good()) {
      const std::string& path = impl_->paths[impl_->current_index];
      uint64_t length = 0;
      s.read(reinterpret_cast<char*>(&length), 8);
      const std::streamsize len_read = s.gcount();
      if (len_read == 8) {
        // Read and verify the length CRC (masked CRC32C over the 8 length
        // bytes) rather than seekg-skipping it.
        uint32_t len_crc = 0;
        s.read(reinterpret_cast<char*>(&len_crc), 4);
        if (s.gcount() != 4) {
          std::fprintf(stderr,
                       "tfrecord: truncated TFRecord (length CRC) in %s at "
                       "offset %lld\n",
                       path.c_str(), static_cast<long long>(offset_));
          return false;
        }
        const uint32_t expected_len_crc =
            MaskedCrc32c(reinterpret_cast<const char*>(&length), sizeof(length));
        if (len_crc != expected_len_crc) {
          std::fprintf(stderr,
                       "tfrecord: length CRC mismatch in %s at offset %lld\n",
                       path.c_str(), static_cast<long long>(offset_));
          return false;
        }

        record_.resize(length);
        s.read(record_.data(), static_cast<std::streamsize>(length));
        if (static_cast<uint64_t>(s.gcount()) != length) {
          // A partial payload read (0 < gcount < length) is genuine
          // truncation, not a clean record boundary: surface it as an
          // error instead of silently advancing to the next shard.
          std::fprintf(stderr,
                       "tfrecord: truncated TFRecord payload in %s at offset "
                       "%lld: read %lld of %llu bytes\n",
                       path.c_str(), static_cast<long long>(offset_),
                       static_cast<long long>(s.gcount()),
                       static_cast<unsigned long long>(length));
          return false;
        } else {
          // Read and verify the payload CRC (masked CRC32C over record_).
          uint32_t data_crc = 0;
          s.read(reinterpret_cast<char*>(&data_crc), 4);
          if (s.gcount() != 4) {
            std::fprintf(stderr,
                         "tfrecord: truncated TFRecord (payload CRC) in %s at "
                         "offset %lld\n",
                         path.c_str(), static_cast<long long>(offset_));
            return false;
          }
          const uint32_t expected_data_crc =
              MaskedCrc32c(record_.data(), record_.size());
          if (data_crc != expected_data_crc) {
            std::fprintf(stderr,
                         "tfrecord: payload CRC mismatch in %s at offset "
                         "%lld\n",
                         path.c_str(), static_cast<long long>(offset_));
            return false;
          }
          offset_ += 8 + 4 + length + 4;
          return true;
        }
      } else if (len_read != 0) {
        // A partial length read at a record boundary is truncation.
        std::fprintf(stderr,
                     "tfrecord: truncated TFRecord (length field) in %s at "
                     "offset %lld: read %lld of 8 bytes\n",
                     path.c_str(), static_cast<long long>(offset_),
                     static_cast<long long>(len_read));
        return false;
      }
      // len_read == 0: clean EOF at a record boundary — fall through to the
      // shard-advance code below.
    }
    // Current shard exhausted (or read failed at boundary). Try next shard.
    if (impl_->current_index + 1 >= impl_->paths.size()) return false;
    impl_->stream.close();
    ++impl_->current_index;
    impl_->stream.clear();
    // pubsetbuf must be called BEFORE open to take effect.
    impl_->stream.rdbuf()->pubsetbuf(impl_->read_buf.data(),
                                     impl_->read_buf.size());
    impl_->stream.open(impl_->paths[impl_->current_index], std::ios::binary);
    if (!impl_->stream.is_open()) return false;
    offset_ = 0;
  }
}

void TFRecordReader::Close() {
  if (impl_) impl_->stream.close();
}

// ---------------------------------------------------------------------------
// TFRecordWriter
// ---------------------------------------------------------------------------
//
// Implementation note (2026-05-01): we used to back this with
// std::ofstream, which buffers writes in a userspace buffer and lets
// the kernel buffer dirty pages indefinitely. On macOS that triggers
// Jetsam after ~137 GB of dirty file-backed memory in a 24h window,
// killing our process mid-WG run. Switched to a raw POSIX fd with
// F_NOCACHE so writes go straight to the disk device without
// accumulating in the kernel page cache. We keep a small userspace
// buffer (kBufBytes) so each fd write is large enough that the SSD
// can actually batch them; no perf regression observed on chr20.

namespace {
constexpr size_t kBufBytes = 1 << 20;  // 1 MiB write coalescing buffer
static_assert(kBufBytes % 4096 == 0,
              "F_NOCACHE requires sector-aligned buffer");
}

struct TFRecordWriter::Impl {
  int fd = -1;
  std::vector<char> buf;
  size_t buf_used = 0;
  bool ok = false;

  explicit Impl(const std::string& path) : buf(kBufBytes) {
    fd = ::open(path.c_str(),
                O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) return;
    // F_NOCACHE: bypass the unified buffer cache. Writes go straight
    // to disk; pages are NOT marked dirty in the kernel's accounting,
    // so Jetsam doesn't accumulate quota. Only available on macOS.
    ::fcntl(fd, F_NOCACHE, 1);
    // Pre-allocate a hint to the FS for sequential write.
    fcntl(fd, F_RDADVISE, 0);  // best-effort; ignored if unsupported
    ok = true;
  }

  ~Impl() {
    FlushBuf();
    if (fd >= 0) ::close(fd);
  }

  bool FlushBuf() {
    if (!ok || buf_used == 0) return ok;
    // BUG FIX (2026-05-10): F_NOCACHE on macOS silently truncates writes
    // that are not multiples of the disk's sector size. Empirically a
    // 155 KiB partial last record at end-of-file got truncated to
    // 139 KiB (= 34 × 4 KiB rounded down) — the kernel writes only the
    // sector-aligned prefix and discards the tail without an error
    // return. This caused 1 record per shard to be lost on close, then
    // the previous TFRecordReader bug (return false on truncated tail)
    // amplified it to 95 % data loss in multi-shard reads.
    //
    // Fix: only the partial-buffer flush (`buf_used < buf.size()`) hits
    // the alignment problem. For full 1-MiB buffer flushes we keep
    // F_NOCACHE on (avoiding macOS Jetsam from dirty-page accounting at
    // WG scale, per the implementation note above). For partial flushes
    // we re-enable the buffered path so the kernel can write any byte
    // count cleanly.
    const bool partial = buf_used < buf.size();
    if (partial && fd >= 0) ::fcntl(fd, F_NOCACHE, 0);

    const char* p = buf.data();
    size_t left = buf_used;
    while (left > 0) {
      ssize_t n = ::write(fd, p, left);
      if (n <= 0) { ok = false; return false; }
      p += n;
      left -= static_cast<size_t>(n);
    }
    buf_used = 0;

    // Re-enable F_NOCACHE for any subsequent full-buffer flushes.
    if (partial && fd >= 0) ::fcntl(fd, F_NOCACHE, 1);
    return true;
  }

  bool Append(const char* data, size_t n) {
    if (!ok) return false;
    while (n > 0) {
      const size_t room = buf.size() - buf_used;
      const size_t take = std::min(n, room);
      std::memcpy(buf.data() + buf_used, data, take);
      buf_used += take;
      data += take;
      n -= take;
      if (buf_used == buf.size()) {
        if (!FlushBuf()) return false;
      }
    }
    return true;
  }
};

TFRecordWriter::TFRecordWriter() = default;
TFRecordWriter::~TFRecordWriter() = default;

std::unique_ptr<TFRecordWriter> TFRecordWriter::New(
    const std::string& path, const std::string& /*compression_type*/) {
  auto impl = std::make_unique<Impl>(path);
  if (!impl->ok) return nullptr;
  auto w = std::unique_ptr<TFRecordWriter>(new TFRecordWriter());
  w->impl_ = std::move(impl);
  return w;
}

bool TFRecordWriter::WriteRecord(std::string_view payload) {
  if (!impl_ || !impl_->ok) return false;
  uint64_t len = payload.size();
  uint32_t len_crc =
      MaskedCrc32c(reinterpret_cast<const char*>(&len), sizeof(len));
  uint32_t data_crc = MaskedCrc32c(payload.data(), len);
  if (!impl_->Append(reinterpret_cast<const char*>(&len), 8)) return false;
  if (!impl_->Append(reinterpret_cast<const char*>(&len_crc), 4)) return false;
  if (!impl_->Append(payload.data(), len)) return false;
  if (!impl_->Append(reinterpret_cast<const char*>(&data_crc), 4)) return false;
  return true;
}

bool TFRecordWriter::WriteRecord(const std::string& payload) {
  return WriteRecord(std::string_view(payload));
}

bool TFRecordWriter::Flush() {
  if (!impl_) return false;
  return impl_->FlushBuf();
}

bool TFRecordWriter::Close() {
  if (!impl_) return true;
  bool ok = impl_->FlushBuf();
  if (impl_->fd >= 0) {
    // Best-effort durability on macOS: flush the device's write cache so the
    // record is on stable storage before we report success.
    ::fcntl(impl_->fd, F_FULLFSYNC, 0);
    if (::close(impl_->fd) != 0) ok = false;
    impl_->fd = -1;
  }
  return ok;
}

}  // namespace deepvariant
