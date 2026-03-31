#ifndef MZ_MATRIX_RR
#define MZ_MATRIX_RR

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "meza/view/base.hpp" // for matrix

#include "quilt/shim.hpp"  // for qs::contains, qs::format
#include "quilt/types.hpp" // for qt::u32, qt::u8, qt::op_t

namespace meza::view::repeated_row
{
using namespace meza::view::base;

template <typename T>
struct repeated_row : public matrix<T> {
	// optional: track which rows and cols have data for quick filtering
	qt::u8 row_has_data_ = 0;
	std::vector<qt::u8> col_has_data_;

	// --------------
	// constructor(s)
	// --------------
	repeated_row(qt::u32 I, qt::u32 J, layout l, T *ptr = nullptr)
	    : matrix<T>(I, J, l, ptr), col_has_data_(J, 0)
	{}

	// ---------
	// modifiers
	// ---------

	void set_has_data(qt::u32 i, qt::u32 j, bool state)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} Index out of range: i={}, j={}",
					   MODULE, i, j));

		row_has_data_ = state;
		col_has_data_[j] = state;
	}

	void set(qt::u32 i, qt::u32 j, T v)
	{
		T *b = this->base;
		std::size_t idx = this->get_idx(i, j);
		b[idx] = v;

		// we never convert a 1 to a 0, so we only need to check if v is
		// non-zero
		if (v != T{})
			this->set_has_data(i, j, true);
	}
};

}; // namespace meza::view::repeated_row

#endif // MZ_MATRIX_RR
