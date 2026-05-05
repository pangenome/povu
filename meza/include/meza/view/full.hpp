#ifndef MZ_MATRIX_FULL
#define MZ_MATRIX_FULL

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "meza/view/base.hpp" // for matrix

#include <log/location.hpp> // for LOG_HERE
#include <quilt/shim.hpp>   // for qs::contains, qs::format
#include <quilt/types.hpp>  // for qt::u32, qt::u8, qt::op_t

namespace meza::view::full
{

using namespace meza::view::base;

template <typename T>
struct full : public matrix<T> {
	std::vector<qt::u8> row_has_data_;
	std::vector<qt::u8> col_has_data_;

	// -----------
	// constructor
	// -----------

	// TODO: should we assume that layout for full is always
	// DenseRowMajor?
	full(qt::u32 I, qt::u32 J, layout l, T *ptr = nullptr)
	    : matrix<T>(I, J, l, ptr), row_has_data_(I, 0), col_has_data_(J, 0)
	{}

	// delete copy constructor to prevent accidental copying of the full
	// matrix
	full(const full &other) = delete;

	// move
	full(full &&other) noexcept
	    : matrix<T>(std::move(other)),
	      row_has_data_(std::move(other.row_has_data_)),
	      col_has_data_(std::move(other.col_has_data_))
	{}

	// -------
	// getters
	// -------

	[[nodiscard]]
	qt::u32 size() const
	{
		return static_cast<qt::u32>(this->rows()) *
		       static_cast<qt::u32>(this->cols());
	}

	[[nodiscard]]
	bool is_row_blank(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(
				qs::format("{} Row index out of range: i={}",
					   LOG_HERE, i));

		return row_has_data_[i] == 0;
	}

	// ---------
	// modifiers
	// ---------

	void comp_metadata()
	{
		const T *e_ptr = this->base;
		const qt::u32 J = this->cols();
		qt::u32 len = this->size();

		for (qt::u32 e_idx{}; e_idx < len; e_idx++, e_ptr++) {
			if (*e_ptr > 0) {
				qt::u32 i = e_idx / J;
				qt::u32 j = e_idx % J;

				this->row_has_data_[i] = 1;
				this->col_has_data_[j] = 1;
			}
		}
	}

	void set_has_data(qt::u32 i, qt::u32 j, bool state)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} Index out of range: i={}, j={}",
					   LOG_HERE, i, j));

		row_has_data_[i] = state;
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

}; // namespace meza::view::full

#endif // MZ_MATRIX_FULL
