/*
 * Copyright 2018 Google LLC.
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
 *
 */

#include "third_party/nucleus/io/reader_base.h"

#include <algorithm>
#include <string>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "third_party/nucleus/core/statusor.h"
#include "tensorflow/core/lib/core/errors.h"

using ::testing::Pointee;
using ::testing::StrEq;

namespace nucleus {

// Forward declaration.
class ToyIterable;

class ToyReader : public Reader {
 private:
  // Underlying container.
  std::vector<StatusOr<string>> toys_;

 public:
  explicit ToyReader(const std::vector<StatusOr<string>>& toys) : toys_(toys) {}

  explicit ToyReader(const std::vector<string>& toys) {
    std::transform(toys.begin(), toys.end(), std::back_inserter(toys_),
                   [](const string& s) { return StatusOr<string>(s); });
  }

  std::shared_ptr<ToyIterable> IterateFrom(int startingPos = 0) {
    if (startingPos < 0) {
      LOG(ERROR) << "ToyIterable: bad starting pos";
      return nullptr;
    } else {
      return MakeIterable<ToyIterable>(this, startingPos);
    }
  }

  friend class ToyIterable;
};

class ToyIterable : public Iterable<string> {
 private:
  uint32 pos_;

 public:
  StatusOr<bool> Next(string* out) override {
    const ToyReader& reader = *static_cast<const ToyReader*>(reader_);
    NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
    if (pos_ < reader.toys_.size()) {
      StatusOr<string> toy_or = reader.toys_[pos_];
      NUCLEUS_RETURN_IF_ERROR(toy_or.status());
      *out = toy_or.ValueOrDie();
      ++pos_;
      return true;
    } else {
      return false;
    }
  }

  ToyIterable(const ToyReader* reader, int startingPos)
      : Iterable(reader), pos_(startingPos) {}

  ~ToyIterable() override {}
};

TEST(ReaderIterableTest, EmptyReaderRange) {
  ToyReader tr0(std::vector<string>{});
  int i = 0;
  for (const StatusOr<string*> toy_status : tr0.IterateFrom(0)) {
    ASSERT_THAT(toy_status, IsOK());
    toy_status.ValueOrDie();
    i++;
  }
  EXPECT_EQ(0, i);
}

TEST(ReaderIterableTest, SupportsBasicIteration) {
  std::vector<string> toys = {"ball", "doll", "house", "legos"};
  std::vector<string> from1 = {"doll", "house", "legos"};
  ToyReader tr(toys);
  auto it = tr.IterateFrom(1);
  auto b = begin(it);
  auto e = end(it);
  std::vector<string> gathered;
  for (auto cur = b; cur != e; ++cur) {
    gathered.push_back(*(*cur).ValueOrDie());
  }
  EXPECT_EQ(from1, gathered);
}

TEST(ReaderIterableTest, SupportsRangeFor) {
  std::vector<string> toys = {"ball", "doll", "house", "legos"};
  std::vector<string> from1 = {"doll", "house", "legos"};
  ToyReader tr(toys);
  std::vector<string> gathered;
  for (const StatusOr<string*> toy : tr.IterateFrom(1)) {
    ASSERT_THAT(toy, IsOK());
    gathered.push_back(*toy.ValueOrDie());
  }
  EXPECT_EQ(from1, gathered);
}

// Ensure that the Iterable Next() interface properly handles an error, for
// example as would be encountered upon parsing a malformed record in a file.
// This interface is used by our Python APIs.
TEST(ReaderIterableTest, IterationHandlesError) {
  ToyReader tr({StatusOr<string>("ball"),
                ::nucleus::Unknown("Malformed record: argybarg"),
                StatusOr<string>("doll")});

  std::shared_ptr<ToyIterable> it = tr.IterateFrom(0);
  StatusOr<bool> not_eof_or;
  string line;

  not_eof_or = it->Next(&line);
  ASSERT_TRUE(not_eof_or.ok() && not_eof_or.ValueOrDie());
  ASSERT_EQ(line, "ball");

  not_eof_or = it->Next(&line);
  ASSERT_THAT(not_eof_or, IsNotOKWithMessage("Malformed record: argybarg"));

  // After initially encountering a failure, successive Next() calls will
  // continue to return the same error--we cannot advance further.
  not_eof_or = it->Next(&line);
  ASSERT_THAT(not_eof_or, IsNotOKWithMessage("Malformed record: argybarg"));
}

// Ensure that C++ iterator interface properly handles an error, for example as
// would be encountered upon parsing a malformed record in a file.
TEST(ReaderIterableTest, CppIterationHandlesError) {
  ToyReader tr({StatusOr<string>("ball"),
                ::nucleus::Unknown("Malformed record: argybarg"),
                StatusOr<string>("doll")});
  auto it = tr.IterateFrom(0);
  auto it_cur = begin(it);
  auto it_end = end(it);

  ASSERT_FALSE(it_cur == it_end);
  ASSERT_THAT(*it_cur, IsOK());
  ASSERT_THAT((*it_cur).ValueOrDie(), Pointee(StrEq("ball")));

  ++it_cur;
  ASSERT_FALSE(it_cur == it_end);
  ASSERT_THAT(*it_cur, IsNotOKWithMessage("Malformed record: argybarg"));

  // We cannot advance any further once an error has been encountered.
  ++it_cur;
  ASSERT_TRUE(it_cur == it_end);
}

TEST(ReaderIterableTest, TestProtectionAgainstMultipleIteration) {
  ToyReader tr({"ball", "doll", "house", "legos"});

  // Scope for RAII auto-destruction of iterable.
  {
    auto it1 = tr.IterateFrom(0);
    auto it2 = tr.IterateFrom(0);

    // The first iterator is good; the second should be null because
    // we detected the attempt to get two concurrent iterators.
    EXPECT_NE(nullptr, it1);
    EXPECT_EQ(nullptr, it2);
  }
  // it1 has died, so we can get a new iterable successfully.
  auto it3 = tr.IterateFrom(0);
  EXPECT_NE(it3, nullptr);

  int i = 0;
  for (const StatusOr<string*> toy : it3) {
    ASSERT_THAT(toy, IsOK());
    (void)toy.ValueOrDie();
    i++;
  }
  EXPECT_EQ(i, 4);
}

TEST(ReaderIterableTest, TestExplicitRelease) {
  ToyReader tr({"ball", "doll", "house", "legos"});
  std::shared_ptr<ToyIterable> it1 = tr.IterateFrom(0);
  EXPECT_NE(it1, nullptr);
  std::shared_ptr<ToyIterable> it2 = tr.IterateFrom(0);
  EXPECT_EQ(it2, nullptr);
  ASSERT_THAT(it1->Release(), IsOK());
  std::shared_ptr<ToyIterable> it3 = tr.IterateFrom(0);
  EXPECT_NE(it3, nullptr);
}

TEST(ReaderIterableTest, TestReaderDiesBeforeIterable) {
  std::shared_ptr<ToyIterable> ti;
  {
    ToyReader tr({"ball", "doll", "house", "legos"});
    ti = tr.IterateFrom(0);
    string s;
    StatusOr<bool> status = ti->Next(&s);
    ASSERT_THAT(status, IsOK());
    EXPECT_TRUE(status.ValueOrDie());
  }
  // tr has been destructed; ti is about to be.  If ti doesn't know
  // the reader is dead, we are likely to crash here.  This can happen
  // in Python since destruction order is non-deterministic.
}

}  // namespace nucleus
