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
 */

#include "third_party/nucleus/testing/protocol-buffer-matchers.h"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/strings/substitute.h"
#include "third_party/nucleus/platform/types.h"
#include "google/protobuf/descriptor.h"
#include "google/protobuf/io/tokenizer.h"
#include "google/protobuf/message.h"
#include "google/protobuf/text_format.h"
#include "google/protobuf/util/message_differencer.h"
#include "re2/re2.h"

namespace nucleus {
namespace internal {

using absl::string_view;
using RegExpStringPiece = absl::string_view;  // copybara

// Utilities.

class StringErrorCollector : public google::protobuf::io::ErrorCollector {
 public:
  explicit StringErrorCollector(string* error_text) : error_text_(error_text) {}

  void AddError(int line, int column, const string& message) override {
    absl::SubstituteAndAppend(error_text_, "$0($1): $2\n",
                              line, column, message.c_str());
  }

  void AddWarning(int line, int column, const string& message) override {
    absl::SubstituteAndAppend(error_text_, "$0($1): $2\n",
                              line, column, message.c_str());
  }

 private:
  string* error_text_;
  StringErrorCollector(const StringErrorCollector&) = delete;
  StringErrorCollector& operator=(const StringErrorCollector&) = delete;
};

bool ParsePartialFromAscii(const string& pb_ascii, google::protobuf::Message* proto,
                           string* error_text) {
  google::protobuf::TextFormat::Parser parser;
  StringErrorCollector collector(error_text);
  parser.RecordErrorsTo(&collector);
  parser.AllowPartialMessage(true);
  return parser.ParseFromString(pb_ascii, proto);
}

// Returns true iff p and q can be compared (i.e. have the same descriptor).
bool ProtoComparable(const google::protobuf::Message& p, const google::protobuf::Message& q) {
  return p.GetDescriptor() == q.GetDescriptor();
}

template <typename Container>
string JoinStringPieces(const Container& strings, string_view separator) {
  std::stringstream stream;
  string_view sep = "";
  for (const string_view str : strings) {
    stream << sep << str;
    sep = separator;
  }
  return stream.str();
}

// Find all the descriptors for the ignore_fields.
std::vector<const google::protobuf::FieldDescriptor*> GetFieldDescriptors(
    const google::protobuf::Descriptor* proto_descriptor,
    const std::vector<string>& ignore_fields) {
  std::vector<const google::protobuf::FieldDescriptor*> ignore_descriptors;
  std::vector<string_view> remaining_descriptors;

  const google::protobuf::DescriptorPool* pool = proto_descriptor->file()->pool();
  for (const string& name : ignore_fields) {
    if (const google::protobuf::FieldDescriptor* field = pool->FindFieldByName(name)) {
      ignore_descriptors.push_back(field);
    } else {
      remaining_descriptors.push_back(name);
    }
  }

  QCHECK(remaining_descriptors.empty())
      << "Could not find fields for proto " << proto_descriptor->full_name()
      << " with fully qualified names: "
      << JoinStringPieces(remaining_descriptors, ",");
  return ignore_descriptors;
}

// Sets the ignored fields corresponding to ignore_fields in differencer. Dies
// if any is invalid.
void SetIgnoredFieldsOrDie(const google::protobuf::Descriptor& root_descriptor,
                           const std::vector<string>& ignore_fields,
                           google::protobuf::util::MessageDifferencer* differencer) {
  if (!ignore_fields.empty()) {
    std::vector<const google::protobuf::FieldDescriptor*> ignore_descriptors =
        GetFieldDescriptors(&root_descriptor, ignore_fields);
    for (std::vector<const google::protobuf::FieldDescriptor*>::iterator it =
             ignore_descriptors.begin();
         it != ignore_descriptors.end(); ++it) {
      differencer->IgnoreField(*it);
    }
  }
}

// A criterion that ignores a field path.
class IgnoreFieldPathCriteria
    : public google::protobuf::util::MessageDifferencer::IgnoreCriteria {
 public:
  explicit IgnoreFieldPathCriteria(
      const std::vector<google::protobuf::util::MessageDifferencer::SpecificField>&
          field_path)
      : ignored_field_path_(field_path) {}

  bool IsIgnored(
      const google::protobuf::Message& message1, const google::protobuf::Message& message2,
      const google::protobuf::FieldDescriptor* field,
      const std::vector<google::protobuf::util::MessageDifferencer::SpecificField>&
          parent_fields) override {
    // The off by one is for the current field.
    if (parent_fields.size() + 1 != ignored_field_path_.size()) {
      return false;
    }
    for (size_t i = 0; i < parent_fields.size(); ++i) {
      const auto& cur_field = parent_fields[i];
      const auto& ignored_field = ignored_field_path_[i];
      // We could compare pointers but it's not guaranteed that descriptors come
      // from the same pool.
      if (cur_field.field->full_name() != ignored_field.field->full_name()) {
        return false;
      }

      // repeated_field[i] is ignored if repeated_field is ignored. To put it
      // differently: if ignored_field specifies an index, we ignore only a
      // field with the same index.
      if (ignored_field.index != -1 && ignored_field.index != cur_field.index) {
        return false;
      }
    }
    return field->full_name() == ignored_field_path_.back().field->full_name();
  }

 private:
  const std::vector<google::protobuf::util::MessageDifferencer::SpecificField>
      ignored_field_path_;
};

namespace {
bool Consume(RegExpStringPiece* s, RegExpStringPiece x) {
  // We use the implementation of ABSL's StartsWith here until we can pick up a
  // dependency on Abseil.
  if (x.empty() ||
      (s->size() >= x.size() && memcmp(s->data(), x.data(), x.size()) == 0)) {
    s->remove_prefix(x.size());
    return true;
  }
  return false;
}
}  // namespace

// Parses a field path and returns individual components.
std::vector<google::protobuf::util::MessageDifferencer::SpecificField>
ParseFieldPathOrDie(const string& relative_field_path,
                    const google::protobuf::Descriptor& root_descriptor) {
  std::vector<google::protobuf::util::MessageDifferencer::SpecificField> field_path;
  // We're parsing a dot-separated list of elements that can be either:
  //   - field names
  //   - extension names
  //   - indexed field names
  // The parser is very permissive as to what is a field name, then we check
  // the field name against the descriptor.

  // Regular parsers. Consume() does not handle optional captures so we split it
  // in two regexps.
  const RE2 field_regex(R"(([^.()[\]]+))");
  const RE2 field_subscript_regex(R"(([^.()[\]]+)\[(\d+)\])");
  const RE2 extension_regex(R"(\(([^)]+)\))");

  RegExpStringPiece input(relative_field_path);
  while (!input.empty()) {
    // Consume a dot, except on the first iteration.
    if (input.size() < relative_field_path.size() && !Consume(&input, ".")) {
      LOG(FATAL) << "Cannot parse field path '" << relative_field_path
                 << "' at offset " << relative_field_path.size() - input.size()
                 << ": expected '.'";
    }
    // Try to consume a field name. If that fails, consume an extension name.
    google::protobuf::util::MessageDifferencer::SpecificField field;
    string name;
    if (RE2::Consume(&input, field_subscript_regex, &name, &field.index) ||
        RE2::Consume(&input, field_regex, &name)) {
      if (field_path.empty()) {
        field.field = root_descriptor.FindFieldByName(name);
        CHECK(field.field) << "No such field '" << name << "' in message '"
                           << root_descriptor.full_name() << "'";
      } else {
        const google::protobuf::util::MessageDifferencer::SpecificField& parent =
            field_path.back();
        field.field = parent.field->message_type()->FindFieldByName(name);
        CHECK(field.field) << "No such field '" << name << "' in '"
                           << parent.field->full_name() << "'";
      }
    } else if (RE2::Consume(&input, extension_regex, &name)) {
      field.field =
          google::protobuf::DescriptorPool::generated_pool()->FindExtensionByName(name);
      CHECK(field.field) << "No such extension '" << name << "'";
      if (field_path.empty()) {
        CHECK(root_descriptor.IsExtensionNumber(field.field->number()))
            << "Extension '" << name << "' does not extend message '"
            << root_descriptor.full_name() << "'";
      } else {
        const google::protobuf::util::MessageDifferencer::SpecificField& parent =
            field_path.back();
        CHECK(parent.field->message_type()->IsExtensionNumber(
            field.field->number()))
            << "Extension '" << name << "' does not extend '"
            << parent.field->full_name() << "'";
      }
    } else {
      LOG(FATAL) << "Cannot parse field path '" << relative_field_path
                 << "' at offset " << relative_field_path.size() - input.size()
                 << ": expected field or extension";
    }
    field_path.push_back(field);
  }

  CHECK(!field_path.empty());
  CHECK(field_path.back().index == -1)
      << "Terminally ignoring fields by index is currently not supported ('"
      << relative_field_path << "')";
  return field_path;
}

// Sets the ignored field paths corresponding to field_paths in differencer.
// Dies if any path is invalid.
void SetIgnoredFieldPathsOrDie(const google::protobuf::Descriptor& root_descriptor,
                               const std::vector<string>& field_paths,
                               google::protobuf::util::MessageDifferencer* differencer) {
  for (const string& field_path : field_paths) {
    differencer->AddIgnoreCriteria(std::make_unique<IgnoreFieldPathCriteria>(
        ParseFieldPathOrDie(field_path, root_descriptor)));
  }
}

// Configures a MessageDifferencer and DefaultFieldComparator to use the logic
// described in comp. The configured differencer is the output of this function,
// but a FieldComparator must be provided to keep ownership clear.
void ConfigureDifferencer(const internal::ProtoComparison& comp,
                          google::protobuf::util::DefaultFieldComparator* comparator,
                          google::protobuf::util::MessageDifferencer* differencer,
                          const google::protobuf::Descriptor* descriptor) {
  differencer->set_message_field_comparison(comp.field_comp);
  differencer->set_scope(comp.scope);
  comparator->set_float_comparison(comp.float_comp);
  comparator->set_treat_nan_as_equal(comp.treating_nan_as_equal);
  differencer->set_repeated_field_comparison(comp.repeated_field_comp);
  SetIgnoredFieldsOrDie(*descriptor, comp.ignore_fields, differencer);
  SetIgnoredFieldPathsOrDie(*descriptor, comp.ignore_field_paths, differencer);
  if (comp.float_comp == internal::kProtoApproximate &&
      (comp.has_custom_margin || comp.has_custom_fraction)) {
    // Two fields will be considered equal if they're within the fraction _or_
    // within the margin. So setting the fraction to 0.0 makes this effectively
    // a "SetMargin". Similarly, setting the margin to 0.0 makes this
    // effectively a "SetFraction".
    comparator->SetDefaultFractionAndMargin(comp.float_fraction,
                                            comp.float_margin);
  }
  differencer->set_field_comparator(comparator);
}

// Returns true iff actual and expected are comparable and match.  The
// comp argument specifies how two are compared.
bool ProtoCompare(const internal::ProtoComparison& comp,
                  const google::protobuf::Message& actual,
                  const google::protobuf::Message& expected) {
  if (!ProtoComparable(actual, expected))
    return false;

  google::protobuf::util::MessageDifferencer differencer;
  google::protobuf::util::DefaultFieldComparator field_comparator;
  ConfigureDifferencer(comp, &field_comparator, &differencer,
                       actual.GetDescriptor());

  // It's important for 'expected' to be the first argument here, as
  // Compare() is not symmetric.  When we do a partial comparison,
  // only fields present in the first argument of Compare() are
  // considered.
  return differencer.Compare(expected, actual);
}

// Describes the types of the expected and the actual protocol buffer.
string DescribeTypes(const google::protobuf::Message& expected,
                     const google::protobuf::Message& actual) {
  return absl::StrCat("whose type should be ",
                      expected.GetDescriptor()->full_name(),
                      " but actually is ", actual.GetDescriptor()->full_name());
}

// Prints the protocol buffer pointed to by proto.
string PrintProtoPointee(const google::protobuf::Message* proto) {
  if (proto == NULL)
    return "";

  return "which points to " + ::testing::PrintToString(*proto);
}

// Describes the differences between the two protocol buffers.
string DescribeDiff(const internal::ProtoComparison& comp,
                    const google::protobuf::Message& actual,
                    const google::protobuf::Message& expected) {
  google::protobuf::util::MessageDifferencer differencer;
  google::protobuf::util::DefaultFieldComparator field_comparator;
  ConfigureDifferencer(comp, &field_comparator, &differencer,
                       actual.GetDescriptor());

  string diff;
  differencer.ReportDifferencesToString(&diff);

  // We must put 'expected' as the first argument here, as Compare()
  // reports the diff in terms of how the protobuf changes from the
  // first argument to the second argument.
  differencer.Compare(expected, actual);

  // Removes the trailing '\n' in the diff to make the output look nicer.
  if (diff.length() > 0 && *(diff.end() - 1) == '\n') {
    diff.erase(diff.end() - 1);
  }

  return "with the difference:\n" + diff;
}

bool ProtoMatcherBase::MatchAndExplain(
    const google::protobuf::Message& arg,
    bool is_matcher_for_pointer,  // true iff this matcher is used to match
                                  // a protobuf pointer.
    ::testing::MatchResultListener* listener) const {
  if (must_be_initialized_ && !arg.IsInitialized()) {
    *listener << "which isn't fully initialized";
    return false;
  }

  const google::protobuf::Message* const expected = CreateExpectedProto(arg, listener);
  if (expected == NULL)
    return false;

  // Protobufs of different types cannot be compared.
  const bool comparable = ProtoComparable(arg, *expected);
  const bool match = comparable && ProtoCompare(comp(), arg, *expected);

  // Explaining the match result is expensive.  We don't want to waste
  // time calculating an explanation if the listener isn't interested.
  if (listener->IsInterested()) {
    const char* sep = "";
    if (is_matcher_for_pointer) {
      *listener << PrintProtoPointee(&arg);
      sep = ",\n";
    }

    if (!comparable) {
      *listener << sep << DescribeTypes(*expected, arg);
    } else if (!match) {
      *listener << sep << DescribeDiff(comp(), arg, *expected);
    }
  }

  DeleteExpectedProto(expected);
  return match;
}

}  // namespace internal
}  // namespace nucleus
