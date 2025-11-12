#ifndef POVU_TYPES_CORE_HPP
#define POVU_TYPES_CORE_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ostream>
#include <sys/types.h>
#include <tuple>

namespace povu::types::core
{
/* type aliases for fixed width types */
using u8 = u_int8_t;
using u32 = u_int32_t;
using status_t = int8_t;			 // return status of a fn
using Time = std::chrono::high_resolution_clock; // C++ timer

// TODO: deprecate and replace id_t and idx_ types with u32
using id_t = u32;
using idx_t = u32;

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

struct slice_t {
private:
	u32 start_;
	u32 len_;

public:
	// -----------
	// constructor
	// -----------
	slice_t(u32 s, u32 l) : start_{s}, len_{l}
	{}

	// -------
	// getter(s)
	// -------
	[[nodiscard]]
	op_t<u32> data() const
	{
		return {start_, len_};
	}

	[[nodiscard]]
	u32 start() const
	{
		return start_;
	}

	[[nodiscard]]
	u32 len() const
	{
		return len_;
	}

	// -------
	// friends
	// -------
	friend bool operator==(const slice_t &s1, const slice_t &s2)
	{
		return s1.start_ == s2.start_ && s1.len_ == s2.len_;
	}

	friend bool operator!=(const slice_t &s1, const slice_t &s2)
	{
		return !(s1 == s2);
	}

	friend std::ostream &operator<<(std::ostream &os, const slice_t &s)
	{
		return os << "(" << s.start_ << "," << s.len_ << ")";
	}
};

} // namespace povu::types::core

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pt = povu::types::core;

#endif
