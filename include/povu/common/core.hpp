#ifndef POVU_TYPES_CORE_HPP
#define POVU_TYPES_CORE_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <sys/types.h>
#include <tuple>

namespace povu::types::core
{
[[deprecated("Deprecated: use u32")]]
typedef u_int32_t id_t;
[[deprecated("Deprecated: use u32")]]
typedef u_int32_t idx_t;

/* type aliases for fixed width types */
using u32 = u_int32_t;
using status_t = int8_t;			 // return status of a fn
using Time = std::chrono::high_resolution_clock; // C++ timer

struct slice_t {
	idx_t start;
	idx_t len;

	slice_t(idx_t start, idx_t len) : start{start}, len{len}
	{}
};

/**
 * an ordered pair type similar to std::pair but with same type on both sides
 * for less typing
 */
template <typename T>
using op_t = std::pair<T, T>;

/**
 * unordered pair with same type on both sides
 * therefore always stores as (min, max)
 * we prefer (a,b) over (l,r) to avoid confusion with (left, right) in unordered
 * pairs
 */
template <typename T>
struct unordered_pair {
	T a_;
	T b_;

	// -----------
	// Constructor
	// -----------
	unordered_pair(T a, T b) : a_(std::min(a, b)), b_(std::max(a, b))
	{}

	// --------------------
	// Comparison operators
	// --------------------

	friend bool operator==(const unordered_pair &up1,
			       const unordered_pair &up2)
	{
		return up1.a_ == up2.a_ && up1.b_ == up2.b_;
	}

	friend bool operator<(const unordered_pair &up1,
			      const unordered_pair &up2)
	{
		return std::tie(up1.a_, up1.b_) < std::tie(up2.a_, up2.b_);
	}

	friend bool operator>(const unordered_pair &up1,
			      const unordered_pair &up2)
	{
		return up2 < up1;
	}
};

template <typename T>
using up_t = unordered_pair<T>;
} // namespace povu::types::core

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pt = povu::types::core;

#endif
