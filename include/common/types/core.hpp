#ifndef POVU_TYPES_CORE_HPP
#define POVU_TYPES_CORE_HPP

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <sys/types.h>
#include <utility>

namespace povu::types {

typedef std::chrono::high_resolution_clock Time; // C++ timer

typedef u_int32_t id_t;
typedef u_int32_t idx_t;
typedef int8_t status_t; // return status of a fn

struct Stride {
  std::size_t start;
  std::size_t length;
};
typedef Stride span;

typedef std::pair<std::size_t, std::size_t> size_t_pair;

/**
 * ordered pair similar to std::pair but with same type on both sides for less typing
 */
template <typename T> struct Pair {
  T first;
  T second;

  friend bool operator<(const Pair &lhs, const Pair &rhs) {
    return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
  }

  // spaceship operator
  friend constexpr auto operator<=>(Pair, Pair) = default;
};

template<typename T>
using op_t = Pair<T>;

/**
 * unordered pair with same type on both sides
 */
template <typename T> struct unordered_pair {
  T l;
  T r;

  unordered_pair(T l, T r) : l(std::min(l, r)), r(std::max(l, r)) {}

  // spaceship operator
  friend constexpr auto operator<=>(const unordered_pair &, const unordered_pair &) = default;
};

template <typename T>
using up_t = unordered_pair<T>;
} // namespace povu::types


#endif
