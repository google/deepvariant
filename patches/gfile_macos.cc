// POSIX replacement for third_party/nucleus/io/gfile.cc.
// Reimplements nucleus::Exists, Glob, ReadableFile, WritableFile
// using std::filesystem + POSIX. No TF runtime.
//
// Implementation note: the original class fields (stream_, file_) are typed
// as our stub types (empty structs). We avoid using them by keeping the real
// implementation state in parallel static maps keyed on `this`.

#include "third_party/nucleus/io/gfile.h"

#include <filesystem>
#include <fstream>
#include <glob.h>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace nucleus {

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

bool Exists(const std::string& filename) {
  return std::filesystem::exists(filename);
}

std::vector<std::string> Glob(const std::string& pattern) {
  std::vector<std::string> results;
  glob_t g{};
  if (::glob(pattern.c_str(), GLOB_TILDE, nullptr, &g) == 0) {
    for (size_t i = 0; i < g.gl_pathc; ++i)
      results.emplace_back(g.gl_pathv[i]);
  }
  globfree(&g);
  return results;
}

// ---------------------------------------------------------------------------
// ReadableFile
// ---------------------------------------------------------------------------

namespace {
struct RFImpl { std::ifstream stream; };
std::mutex rf_mu;
std::unordered_map<ReadableFile*, std::unique_ptr<RFImpl>> rf_map;
}  // namespace

ReadableFile::ReadableFile()  = default;
ReadableFile::~ReadableFile() {
  std::lock_guard<std::mutex> lk(rf_mu);
  rf_map.erase(this);
}

std::unique_ptr<ReadableFile> ReadableFile::New(const std::string& filename) {
  auto impl = std::make_unique<RFImpl>();
  impl->stream.open(filename);
  if (!impl->stream.is_open()) return nullptr;
  auto f = std::unique_ptr<ReadableFile>(new ReadableFile());
  {
    std::lock_guard<std::mutex> lk(rf_mu);
    rf_map[f.get()] = std::move(impl);
  }
  return f;
}

bool ReadableFile::Readline(std::string* s) {
  std::lock_guard<std::mutex> lk(rf_mu);
  auto it = rf_map.find(this);
  if (it == rf_map.end()) return false;
  return static_cast<bool>(std::getline(it->second->stream, *s));
}

void ReadableFile::Close() {
  std::lock_guard<std::mutex> lk(rf_mu);
  auto it = rf_map.find(this);
  if (it != rf_map.end()) it->second->stream.close();
}

// ---------------------------------------------------------------------------
// WritableFile
// ---------------------------------------------------------------------------

namespace {
struct WFImpl { std::ofstream stream; };
std::mutex wf_mu;
std::unordered_map<WritableFile*, std::unique_ptr<WFImpl>> wf_map;
}  // namespace

WritableFile::WritableFile()  = default;
WritableFile::~WritableFile() {
  std::lock_guard<std::mutex> lk(wf_mu);
  wf_map.erase(this);
}

std::unique_ptr<WritableFile> WritableFile::New(const std::string& filename) {
  auto impl = std::make_unique<WFImpl>();
  impl->stream.open(filename);
  if (!impl->stream.is_open()) return nullptr;
  auto f = std::unique_ptr<WritableFile>(new WritableFile());
  {
    std::lock_guard<std::mutex> lk(wf_mu);
    wf_map[f.get()] = std::move(impl);
  }
  return f;
}

bool WritableFile::Write(const std::string& s) {
  std::lock_guard<std::mutex> lk(wf_mu);
  auto it = wf_map.find(this);
  if (it == wf_map.end()) return false;
  it->second->stream.write(s.data(), static_cast<std::streamsize>(s.size()));
  return it->second->stream.good();
}

void WritableFile::Close() {
  std::lock_guard<std::mutex> lk(wf_mu);
  auto it = wf_map.find(this);
  if (it != wf_map.end()) it->second->stream.close();
}

}  // namespace nucleus
