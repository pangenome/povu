#ifndef POVU_TYPES_CORE_HPP
#define POVU_TYPES_CORE_HPP

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <sys/types.h>
#include <utility>
#include <tuple>
#include <algorithm>


#include <fmt/format.h>
#include <fmt/color.h>

#include "./compat.hpp"

namespace povu::types {

// =======
// logging
// =======

#define FN() pv_cmp::format("[{}::{}]", MODULE, __func__)

#define DEBUG_PRINT(format, ...)                                               \
  fprintf(stderr, "%s:%d: " format, __FILE__, __LINE__, ##__VA_ARGS__)

/**
 * Generic logging macro: prints a label (e.g. "ERR") with color and formatted
 * output.
 */
#define LOG(label, color, fmt_str, ...)                                        \
  do {                                                                         \
    fmt::print(stderr, fmt::fg(color) | fmt::emphasis::bold, "{} {} ", label,  \
               FN());                                                          \
    fmt::print(stderr, fmt_str, ##__VA_ARGS__);                                \
    fmt::print(stderr, "\n");                                                  \
  } while (false)

#define ERR(fmt_str, ...)                                                      \
  LOG("ERR", fmt::color::crimson, fmt_str, ##__VA_ARGS__)

#define WARN(fmt_str, ...)                                                     \
  LOG("WARN", fmt::color::yellow, fmt_str, ##__VA_ARGS__)


// =====
// types
// =====

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


#endif
