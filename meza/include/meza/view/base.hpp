#ifndef MZ_MATRIX_VIEW_BASE
#define MZ_MATRIX_VIEW_BASE

#include <cstddef>
#include <stdexcept>
#include <string_view>

#include "meza/shared/shared.hpp" // for layout

#include "quilt/shim.hpp"  // for qs::contains, qs::format
#include "quilt/types.hpp" // for qt::u32, qt::u8, qt::op_t

namespace meza::view::base
{
inline constexpr std::string_view MODULE = "meza::matrix::non_own";

using layout = meza::shared::layout;

template <typename T>
struct matrix {
	layout lo;
	qt::u32 I_;
	qt::u32 J_;
	T *base = nullptr;

	matrix(qt::u32 I, qt::u32 J, layout l, T *ptr = nullptr)
	    : lo(l), I_(I), J_(J), base(ptr)
	{}

	// -------
	// getters
	// --------

	[[nodiscard]]
	qt::u32 rows() const
	{
		return I_;
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return J_;
	}

	template <typename V>
	void copy_slice(const V *src_ptr, qt::u32 len, qt::u32 repeat)
	{
		// pointers to mutate
		const V *s; //
		T *t = this->base;

		for (qt::u32 r{}; r < repeat; r++) {
			s = src_ptr; // reset to the src start
			for (qt::u32 i{}; i < len; i++, s++, t++)
				*t = static_cast<T>(*s);
		}
	}

	const T *data() const
	{
		return this->base;
	}

	[[nodiscard]]
	const T &at(qt::u32 i, qt::u32 j) const
	{
		std::size_t idx = this->get_idx(i, j);
		return this->base[idx];
	}

	[[nodiscard]]
	std::size_t get_idx(qt::u32 i, qt::u32 j) const
	{
		qt::u32 stride = this->cols();
		if (this->lo == layout::DenseRowMajor) {
			return (size_t)i * stride + j;
		}
		else if (this->lo == layout::LowerSymmetricSquare) {
			if (i < j)
				throw std::out_of_range(
					qs::format("{} Index out of range for "
						   "LowerSymmetricSquare: "
						   "i={}, j={}",
						   MODULE, i, j));
			return (size_t)i * (i + 1) / 2 + j;
		}
		else if (this->lo == layout::RepeatedRow) {
			return j; // all rows are the same, so just return the
				  // column index
		}
		else {
			throw std::logic_error(
				qs::format("{} Invalid layout: {}", MODULE,
					   static_cast<int>(lo)));
		}
	}
};

} // namespace meza::view::base
#endif // MZ_MATRIX_VIEW_BASE
