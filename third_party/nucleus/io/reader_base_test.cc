/*
 * Copyright 2018 Google Inc.
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

#include <iostream>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/vendor/status_matchers.h"

using tensorflow::string;

namespace nucleus {

// Forward declaration.
class ToyIterable;

class ToyReader : public Reader {
 private:
  // Underlying container.
  std::vector<string> toys_;

 public:
  explicit ToyReader(const std::vector<string>& toys)
      : toys_(toys)
  {}

  std::shared_ptr<ToyIterable> IterateFrom(int startingPos = 0)  {
    if (startingPos < 0) {
      LOG(ERROR) << "ToyIterable: bad starting pos";
      return nullptr;
    } else {
      return MakeIterable<ToyIterable>(this, toys_, startingPos);
    }
  }
};

class ToyIterable : public Iterable<string>  {
 private:
  const std::vector<string>& toys_;
  uint32_t pos_;

 public:
  StatusOr<bool> Next(string* out) override {
    TF_RETURN_IF_ERROR(CheckIsAlive());
    if (pos_ < toys_.size()) {
      *out = toys_[pos_];
      ++pos_;
      return true;
    } else {
      return false;
    }
  }

  ToyIterable(const ToyReader* reader,
              const std::vector<string>& toys,
              int startingPos)
      : Iterable(reader),
        toys_(toys),
        pos_(startingPos)
  {}

  ~ToyIterable() override {}
};



TEST(ReaderIterableTest, EmptyReaderRange) {
  ToyReader tr0({});
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
