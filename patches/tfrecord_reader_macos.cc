// POSIX replacement for third_party/nucleus/io/tfrecord_reader.cc.
// TFRecord format: [uint64_le length][uint32_le masked_crc32c(len)]
//                  [bytes payload][uint32_le masked_crc32c(payload)]
// Uncompressed only. Masked CRC32Cs are verified against the data.

#include "third_party/nucleus/io/tfrecord_reader.h"

#include <cstdint>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>

#include "absl/crc/crc32c.h"
#include "absl/log/log.h"

namespace nucleus {

namespace {
// Mask delta and helper must match those used by the writer
// (tfrecord_writer_macos.cc) so written CRCs verify here.
constexpr uint32_t kMaskDelta = 0xa282ead8UL;
uint32_t MaskedCrc32c(const char* data, size_t n) {
  uint32_t crc = static_cast<uint32_t>(
      absl::ComputeCrc32c(std::string_view(data, n)));
  return ((crc >> 15) | (crc << 17)) + kMaskDelta;
}

struct TFRRImpl {
  std::ifstream stream;
  std::string path;
  explicit TFRRImpl(const std::string& p)
      : stream(p, std::ios::binary), path(p) {}
};
std::mutex mu;
std::unordered_map<TFRecordReader*, std::unique_ptr<TFRRImpl>> impls;
}  // namespace

TFRecordReader::TFRecordReader() : offset_(0) {}
TFRecordReader::~TFRecordReader() {
  std::lock_guard<std::mutex> lk(mu);
  impls.erase(this);
}

// static
std::unique_ptr<TFRecordReader> TFRecordReader::New(
    const std::string& filename, const std::string& /*compression_type*/) {
  auto impl = std::make_unique<TFRRImpl>(filename);
  if (!impl->stream.is_open()) return nullptr;
  auto r = std::unique_ptr<TFRecordReader>(new TFRecordReader());
  {
    std::lock_guard<std::mutex> lk(mu);
    impls[r.get()] = std::move(impl);
  }
  return r;
}

bool TFRecordReader::GetNext() {
  std::lock_guard<std::mutex> lk(mu);
  auto it = impls.find(this);
  if (it == impls.end()) return false;
  auto& s = it->second->stream;
  const std::string& path = it->second->path;
  if (!s.good()) return false;

  // Read the 8-byte little-endian length. A zero-byte read at a record
  // boundary is a clean EOF; a partial read (1..7 bytes) is truncation.
  uint64_t length = 0;
  s.read(reinterpret_cast<char*>(&length), 8);
  std::streamsize length_read = s.gcount();
  if (length_read == 0) return false;  // clean EOF
  if (length_read != 8) {
    LOG(ERROR) << "Truncated TFRecord (short length field, " << length_read
               << "/8 bytes) in " << path;
    return false;
  }

  // Read and verify the masked CRC32C of the 8-byte length field.
  uint32_t length_crc = 0;
  s.read(reinterpret_cast<char*>(&length_crc), 4);
  if (s.gcount() != 4) {
    LOG(ERROR) << "Truncated TFRecord (short length CRC) in " << path;
    return false;
  }
  uint32_t expected_length_crc =
      MaskedCrc32c(reinterpret_cast<const char*>(&length), 8);
  if (length_crc != expected_length_crc) {
    LOG(ERROR) << "Corrupt TFRecord (length CRC mismatch) in " << path;
    return false;
  }

  // Read the payload. A short payload read is truncation, never clean EOF
  // (we are mid-record after consuming a valid length field).
  record_.resize(length);
  s.read(record_.data(), static_cast<std::streamsize>(length));
  if (static_cast<uint64_t>(s.gcount()) != length) {
    LOG(ERROR) << "Truncated TFRecord (short payload, " << s.gcount() << "/"
               << length << " bytes) in " << path;
    return false;
  }

  // Read and verify the masked CRC32C of the payload.
  uint32_t data_crc = 0;
  s.read(reinterpret_cast<char*>(&data_crc), 4);
  if (s.gcount() != 4) {
    LOG(ERROR) << "Truncated TFRecord (short payload CRC) in " << path;
    return false;
  }
  uint32_t expected_data_crc = MaskedCrc32c(record_.data(), length);
  if (data_crc != expected_data_crc) {
    LOG(ERROR) << "Corrupt TFRecord (payload CRC mismatch) in " << path;
    return false;
  }

  offset_ += 8 + 4 + length + 4;
  return true;
}

void TFRecordReader::Close() {
  std::lock_guard<std::mutex> lk(mu);
  auto it = impls.find(this);
  if (it != impls.end()) it->second->stream.close();
}

}  // namespace nucleus
