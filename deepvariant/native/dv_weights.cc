#include "deepvariant/native/dv_weights.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstring>
#include <utility>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr char kMagic[4] = {'D', 'V', 'W', '1'};

// Read N little-endian bytes at `data` interpreted as the native integer
// of the same size.  M-series is little-endian so memcpy is fine.
template <typename T>
inline T ReadLE(const uint8_t* p) {
  T v;
  std::memcpy(&v, p, sizeof(T));
  return v;
}

}  // namespace

DvwWeights::DvwWeights() = default;

DvwWeights::~DvwWeights() {
  if (map_addr_ != nullptr && map_size_ > 0) {
    munmap(map_addr_, map_size_);
  }
}

std::unique_ptr<DvwWeights> DvwWeights::Open(const std::string& path) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    LOG(ERROR) << "DvwWeights: open(" << path << ") failed";
    return nullptr;
  }
  struct stat st;
  if (::fstat(fd, &st) != 0 || st.st_size < 12) {
    LOG(ERROR) << "DvwWeights: fstat(" << path << ") failed or file too small";
    ::close(fd);
    return nullptr;
  }
  void* addr = ::mmap(nullptr, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
  ::close(fd);
  if (addr == MAP_FAILED) {
    LOG(ERROR) << "DvwWeights: mmap(" << path << ") failed";
    return nullptr;
  }

  auto* base = static_cast<const uint8_t*>(addr);
  // Header.
  if (std::memcmp(base, kMagic, 4) != 0) {
    LOG(ERROR) << "DvwWeights: bad magic in " << path;
    munmap(addr, st.st_size);
    return nullptr;
  }
  const uint32_t version = ReadLE<uint32_t>(base + 4);
  const uint32_t n_tensors = ReadLE<uint32_t>(base + 8);
  if (version != 1u) {
    LOG(ERROR) << "DvwWeights: unsupported version " << version;
    munmap(addr, st.st_size);
    return nullptr;
  }

  auto out = std::unique_ptr<DvwWeights>(new DvwWeights());
  out->map_addr_ = addr;
  out->map_size_ = static_cast<size_t>(st.st_size);
  out->version_ = version;
  out->names_.reserve(n_tensors);
  out->by_name_.reserve(n_tensors);

  // Walk the per-tensor table to find the payload start.
  size_t p = 12;
  // First pass to compute table size, then parse entries with payload base.
  size_t table_start = p;
  for (uint32_t t = 0; t < n_tensors; ++t) {
    if (p + 4 > out->map_size_) goto truncated;
    const uint32_t name_len = ReadLE<uint32_t>(base + p);
    p += 4;
    if (p + name_len + 2 > out->map_size_) goto truncated;
    p += name_len;             // name bytes
    p += 1;                    // dtype
    const uint8_t ndim = base[p];
    p += 1;
    if (p + 4u * ndim + 16 > out->map_size_) goto truncated;
    p += 4u * ndim;            // shape
    p += 16;                   // offset + n_bytes
  }
  // p now points at the start of the payload.
  {
    const size_t payload_base = p;

    // Second pass: actually populate by_name_.
    p = table_start;
    for (uint32_t t = 0; t < n_tensors; ++t) {
      const uint32_t name_len = ReadLE<uint32_t>(base + p);
      p += 4;
      std::string name(reinterpret_cast<const char*>(base + p), name_len);
      p += name_len;
      const uint8_t dtype = base[p];
      p += 1;
      const uint8_t ndim = base[p];
      p += 1;
      std::vector<uint32_t> shape(ndim);
      for (uint8_t d = 0; d < ndim; ++d) {
        shape[d] = ReadLE<uint32_t>(base + p);
        p += 4;
      }
      const uint64_t offset = ReadLE<uint64_t>(base + p);
      p += 8;
      const uint64_t n_bytes = ReadLE<uint64_t>(base + p);
      p += 8;
      if (dtype != 1u) {
        LOG(ERROR) << "DvwWeights: tensor " << name
                   << " has unsupported dtype " << static_cast<int>(dtype);
        munmap(addr, st.st_size);
        return nullptr;
      }
      const size_t abs_offset = payload_base + offset;
      if (abs_offset + n_bytes > out->map_size_) {
        LOG(ERROR) << "DvwWeights: tensor " << name
                   << " spills past end of file";
        munmap(addr, st.st_size);
        return nullptr;
      }
      DvwTensor tensor;
      tensor.data = reinterpret_cast<const float*>(base + abs_offset);
      tensor.shape = std::move(shape);
      tensor.n_bytes = n_bytes;
      tensor.n_elements = n_bytes / sizeof(float);
      out->names_.push_back(name);
      out->by_name_.emplace(std::move(name), std::move(tensor));
    }
  }
  return out;

truncated:
  LOG(ERROR) << "DvwWeights: truncated file " << path;
  munmap(addr, st.st_size);
  return nullptr;
}

const DvwTensor* DvwWeights::Get(const std::string& name) const {
  auto it = by_name_.find(name);
  return it == by_name_.end() ? nullptr : &it->second;
}

}  // namespace deepvariant
