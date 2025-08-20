#ifndef POVU_TYPES_CORE_HPP
#define POVU_TYPES_CORE_HPP

#include <chrono>
#include <cstdint>
#include <sys/types.h>
#include <tuple>
#include <algorithm>

namespace povu::types::core {

typedef std::chrono::high_resolution_clock Time; // C++ timer

typedef u_int32_t id_t;
typedef u_int32_t idx_t;
typedef int8_t status_t; // return status of a fn

/**
 * ordered pair similar to std::pair but with same type on both sides for less typing
 */
template <typename T> struct Pair {
  T first;
  T second;

  // --------------------
  // Comparison operators
  // --------------------

  friend bool operator<(const Pair &lhs, const Pair &rhs) {
    return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
  }

  friend bool operator==(const Pair &lhs, const Pair &rhs) {
    return std::tie(lhs.first, lhs.second) == std::tie(rhs.first, rhs.second);
  }

  friend bool operator>(const Pair &lhs, const Pair &rhs) {
    return rhs < lhs;
  }
};

template<typename T>
using op_t = Pair<T>;

/**
 * unordered pair with same type on both sides
 * therefore always stores as (min, max)
 * we prefer (a,b) over (l,r) to avoid confusion with (left, right) in unordered pairs
 */
template <typename T> struct unordered_pair {
  T a_;
  T b_;

  // -----------
  // Constructor
  // -----------
  unordered_pair(T a, T b) : a_(std::min(a, b)), b_(std::max(a, b)) {}

  // --------------------
  // Comparison operators
  // --------------------

  friend bool operator==(const unordered_pair &up1, const unordered_pair &up2) {
    return up1.a_ == up2.a_ && up1.b_ == up2.b_;
  }

  friend bool operator<(const unordered_pair &up1, const unordered_pair &up2) {
    return std::tie(up1.a_, up1.b_) < std::tie(up2.a_, up2.b_);
  }

  friend bool operator>(const unordered_pair &up1, const unordered_pair &up2) {
    return up2 < up1;
  }
};

template <typename T>
using up_t = unordered_pair<T>;
} // namespace povu::types

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pt = povu::types::core;

#endif
