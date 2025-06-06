/*
 * Copyright 2025 Google LLC.
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
#ifndef LEARNING_GENOMICS_DEEPVARIANT_DISTRIBUTION_FUNCTOR_H_
#define LEARNING_GENOMICS_DEEPVARIANT_DISTRIBUTION_FUNCTOR_H_

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <numeric>
#include <optional>
#include <stack>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"

namespace learning::genomics::deepvariant::distribution_functor {

/// A distribution over a finite set of values.
// The main use of this class is to define an input distribution, map it through
// a function to compute an output distribution.
template <typename T>
class Distribution {
 public:
  using value_type = T;

  using WeightMap = absl::flat_hash_map<value_type, uint64_t>;

  const WeightMap& get_weight_map() const { return weight_map_; }

  const uint64_t& get_total_weight() const { return total_weight_; }

  static auto FromWeightMap(const WeightMap& weight_map) {
    typename Distribution<T>::WeightMap weight_map_processed;
    uint64_t total_weight = 0;
    for (auto& [t, p] : weight_map) {
      if (p > 0) {
        weight_map_processed[t] = p;
        total_weight += p;
      }
    }
    return Distribution(weight_map_processed, total_weight);
  }

  friend bool operator==(const Distribution& a, const Distribution& b) {
    return a.weight_map_ == b.weight_map_;
  }

  friend bool operator!=(const Distribution& a, const Distribution& b) {
    return !(a == b);
  }

 private:
  explicit Distribution(const WeightMap& weight_map, uint64_t total_weight)
      : weight_map_(weight_map), total_weight_(total_weight) {
    normalize();
  }

  void normalize() {
    uint64_t factor = total_weight_;
    for (auto& [_, p] : weight_map_) {
      factor = std::gcd(factor, p);
    }
    for (auto& [_, p] : weight_map_) {
      p /= factor;
    }
    total_weight_ /= factor;
  }

  WeightMap weight_map_;
  uint64_t total_weight_;
};

/// Represents a (potentially infinite) product of distributions over the
// type T with parameters of type U. We can still map this (potentially
// infinite) product measure through a function to get a finite distribution,
// since this function can use only finite randomness.
template <typename T, typename... U>
class DistributionGenerator {
 public:
  explicit DistributionGenerator(std::function<Distribution<T>(U...)> g)
      : g_(g) {}
  Distribution<T> operator()(U... u) { return g_(u...); }

 private:
  std::function<Distribution<T>(U...)> g_;
};

/// Traces the calls to a DistributionGenerator, and computes the resulting
// output distribution.
template <typename T, typename R, typename... U>
class DistributionTracer {
 public:
  DistributionTracer(std::function<R(std::function<T(U...)>)> f,
                     DistributionGenerator<T, U...> dist_gen)
      : f_(f), dist_gen_(dist_gen) {}

  Distribution<R> Trace() {
    nodes_to_explore_ = std::stack<std::shared_ptr<TracerNode>>();
    replay_stack_ = std::stack<T>();
    current_node_ = nullptr;
    std::function<T(U...)> callable_provider = [&](U... u) {
      return this->provide_next_response(std::forward<U>(u)...);
    };

    absl::flat_hash_map<R, std::pair<uint64_t, uint64_t>> responses;
    R r = f_(callable_provider);
    // If f_ doesn't actually use the provider, the response is deterministic.
    if (current_node_ == nullptr) {
      return Distribution<R>::FromWeightMap({{r, 1}});
    }

    // At this point, we know f_ uses the provider, so we can start exploring.
    responses[r] = {current_node_->num_, current_node_->denom_};

    while (!nodes_to_explore_.empty()) {
      current_node_ = nodes_to_explore_.top();
      if (current_node_->domain_it_ ==
          current_node_->distribution_->get_weight_map().cend()) {
        nodes_to_explore_.pop();
        continue;
      }
      std::shared_ptr<TracerNode> helper_node = current_node_;
      while (helper_node->response_.has_value()) {
        replay_stack_.push(helper_node->response_.value());
        helper_node = helper_node->parent_;
      }
      R r = f_(callable_provider);
      auto response_it = responses.find(r);
      if (response_it == responses.end()) {
        responses[r] = {current_node_->num_, current_node_->denom_};
      } else {
        auto [prev_p, prev_q] = response_it->second;
        uint64_t numerator =
            current_node_->num_ * prev_q + prev_p * current_node_->denom_;
        uint64_t denominator = current_node_->denom_ * prev_q;
        uint64_t common = std::gcd(numerator, denominator);
        responses[r] = {numerator / common, denominator / common};
      }
    }
    // Compute a common denominator for all the responses.
    uint64_t total_weight = 1;
    for (auto& [_, p_q] : responses) {
      auto common = std::gcd(p_q.second, total_weight);
      total_weight = common * (p_q.second / common) * (total_weight / common);
    }

    typename Distribution<R>::WeightMap resulting_weight_map;
    for (auto& [r, p_q] : responses) {
      resulting_weight_map[r] = p_q.first * (total_weight / p_q.second);
    }
    return Distribution<R>::FromWeightMap(resulting_weight_map);
  }

  T provide_next_response(U... u) {
    if (!replay_stack_.empty()) {
      T t = replay_stack_.top();
      replay_stack_.pop();
      return t;
    }

    if (current_node_ == nullptr) {
      // Create a root node
      current_node_ = std::make_shared<TracerNode>();
    }
    // Set query if not already set, and push the node to the explore stack.
    if (!current_node_->distribution_.has_value()) {
      current_node_->SetDistribution(std::move(dist_gen_(u...)));
      nodes_to_explore_.push(current_node_);
    }

    // Find a suitable response.
    auto [t, w] = *(current_node_->domain_it_++);

    auto gcd1 = std::gcd(current_node_->num_,
                         current_node_->distribution_->get_total_weight());
    auto gcd2 = std::gcd(w, current_node_->denom_);

    uint64_t p = (current_node_->num_ / gcd1) * (w / gcd2);
    uint64_t q = (current_node_->denom_ / gcd2) *
                 (current_node_->distribution_->get_total_weight() / gcd1);

    // Create a child node and update the current node.
    current_node_ = std::make_shared<TracerNode>(current_node_, t, p, q);
    return t;
  }

 private:
  class TracerNode {
   public:
    explicit TracerNode(std::shared_ptr<TracerNode> parent = nullptr,
                        std::optional<T> response = std::nullopt,
                        uint64_t num = 1, uint64_t denom = 1)
        : parent_(parent), response_(response), num_(num), denom_(denom) {}

    void SetDistribution(Distribution<T> dist) {
      distribution_ = dist;
      domain_it_ = distribution_->get_weight_map().cbegin();
    }

    std::optional<Distribution<T>> distribution_;
    typename Distribution<T>::WeightMap::const_iterator domain_it_;
    std::shared_ptr<TracerNode> parent_;
    std::optional<T> response_;
    uint64_t num_, denom_;
  };
  const std::function<R(std::function<T(U...)>)> f_;
  const std::function<Distribution<T>(U...)> dist_gen_;
  std::stack<std::shared_ptr<TracerNode>> nodes_to_explore_;
  std::shared_ptr<TracerNode> current_node_ = nullptr;
  std::stack<T> replay_stack_;
};

//// Simple constructors for distributions.

// Returns a distribution over a single value.
template <typename T>
Distribution<T> unit(T t) {
  return Distribution<T>::FromWeightMap({{t, 1}});
}

// Returns a distribution over the given domain.
template <typename T>
Distribution<T> uniform(std::vector<T> domain) {
  typename Distribution<T>::WeightMap new_w;
  for (const auto& u : domain) {
    new_w[u] = 1;
  }
  return Distribution<T>::FromWeightMap(new_w);
}

//// Operators to transform distributions.

/// Given a distribution over a type T, and a function f that takes a value of
// type T and returns a distribution over a new type U, returns the overall
// output distribution over the type U.
template <typename F, typename T>
auto dist_bind(const Distribution<T>& dist, F f)
    -> Distribution<typename std::invoke_result_t<F, T>::value_type> {
  using OutputDistributionType = std::invoke_result_t<F, T>;
  static_assert(std::is_constructible_v<
                    Distribution<typename OutputDistributionType::value_type>,
                    OutputDistributionType>,
                "f must return a type convertible to Distribution<U>.");
  absl::flat_hash_map<typename OutputDistributionType::value_type,
                      std::pair<uint64_t, uint64_t>>
      new_probs;
  // F is a function from T to Distribution<U>.
  for (const auto& [t, p] : dist.get_weight_map()) {
    auto local_dist = f(t);
    typename OutputDistributionType::WeightMap intermediate_w =
        local_dist.get_weight_map();
    for (const auto& [u, q] : intermediate_w) {
      auto a = std::gcd(p, local_dist.get_total_weight());
      auto b = std::gcd(q, dist.get_total_weight());
      auto product_of_probs = std::make_pair(
          (p / a) * (q / b),
          (dist.get_total_weight() / b) * (local_dist.get_total_weight() / a));
      auto new_probs_it = new_probs.find(u);
      if (new_probs_it == new_probs.end()) {
        new_probs[u] = product_of_probs;
      } else {
        auto [p, q] = product_of_probs;
        auto [prev_p, prev_q] = new_probs_it->second;
        new_probs_it->second = {prev_p * q + p * prev_q, prev_q * q};
      }
    }
  }

  // find a common denominator for all the probabilities.
  uint64_t total_weight = 1;
  for (auto& [_, p] : new_probs) {
    auto common = std::gcd(p.second, total_weight);
    total_weight = p.second * (total_weight / common);
  }

  typename OutputDistributionType::WeightMap new_w;
  for (auto& [u, p] : new_probs) {
    new_w[u] = p.first * (total_weight / p.second);
  }

  return OutputDistributionType::FromWeightMap(new_w);
}

// Given a distribution over a type T, and a function f that takes a value of
// type T and returns a new value of type U, returns the distribution over the
// new type U.
template <typename F, typename T>
auto dist_map(const Distribution<T>& dist, F f)
    -> Distribution<typename std::invoke_result_t<F, T>> {
  using OutputType = std::invoke_result_t<F, T>;
  return dist_bind(dist,
                   [f](T t) -> Distribution<OutputType> { return unit(f(t)); });
}

// Given a distribution generator over a type T with parameters of type U...,
// and a function f that takes a sample generator of type T(U...) and returns a
// value of type R, returns the output distribution over the type R.
template <typename F, typename T, typename... U>
auto dist_map(const DistributionGenerator<T, U...>& dist, const F& f)
    -> Distribution<typename std::invoke_result_t<F, std::function<T(U...)>>> {
  using OutputType = std::invoke_result_t<F, std::function<T(U...)>>;
  return DistributionTracer<T, OutputType, U...>(f, dist).Trace();
}

}  // namespace learning::genomics::deepvariant::distribution_functor

#endif  // LEARNING_GENOMICS_DEEPVARIANT_DISTRIBUTION_FUNCTOR_H_
