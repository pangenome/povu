#ifndef MZ_MATRIX_VIEW
#define MZ_MATRIX_VIEW

#include <iostream>
#include <stdexcept>
#include <vector>

#include <log/location.hpp> // for LOG_HERE
#include <quilt/shim.hpp>   // for qs::contains, qs::format
#include <quilt/types.hpp>  // for qt::u32, qt::u8, qt::op_t

#include "meza/shared/shared.hpp"     // for layout
#include "meza/view/full.hpp"	      // for full
#include "meza/view/repeated_row.hpp" // for repeated_row

namespace meza::view
{

// -----
// type aliases
// -----

template <typename T>
using full_matrix = meza::view::full::full<T>;

template <typename T>
using rr_matrix = meza::view::repeated_row::repeated_row<T>;

// -----
// metadata structs
// -----

template <typename U>
struct col_names {
	std::vector<U> names_;
};

template <typename W>
struct row_names {
	std::vector<W> names_;
};

template <typename Derived, typename MatrixType, typename T, typename U,
	  typename W>
struct matrix_wrapper {
	MatrixType base_matrix_;
	col_names<U> col_names_;
	row_names<W> row_names_;
	bool is_tangled_ = false;
	qt::u32 max_depth_ = 0;

	// --------------
	// constructor(s)
	// --------------

	// Add this constructor
	explicit matrix_wrapper(MatrixType &&m) : base_matrix_(std::move(m))
	{}

	// ----------
	// forwarding
	// ----------
	const MatrixType &base() const
	{
		return base_matrix_;
	}

	MatrixType &base_mut()
	{
		return base_matrix_;
	}

	template <typename V>
	void copy_slice(const V *src_ptr, qt::u32 len, qt::u32 repeat = 1)
	{
		base_mut().copy_slice(src_ptr, len, repeat);
	}

	void comp_meta_sync()
	{
		base_mut().comp_metadata();
	}

	void set_value(qt::u32 i, qt::u32 j, T value)
	{
		base_mut().set(i, j, value);
	}

	[[nodiscard]] T get_value(qt::u32 i, qt::u32 j) const
	{
		return base().at(i, j);
	}

	[[nodiscard]] qt::u32 rows() const
	{
		return base().rows();
	}

	[[nodiscard]] qt::u32 cols() const
	{
		return base().cols();
	}

	// -------
	// getters
	// -------

	const std::vector<U> &get_col_names() const
	{
		return col_names_.names_;
	}

	const std::vector<W> &get_row_names() const
	{
		return row_names_.names_;
	}

	[[nodiscard]]
	qt::u32 get_max_depth() const
	{
		return max_depth_;
	}

	[[nodiscard]]
	bool is_tangled() const
	{
		return is_tangled_;
	}

	// -------
	// setters
	// -------

	void set_max_depth(qt::u32 d)
	{
		max_depth_ = d;
	}

	void set_tangled(bool v)
	{
		is_tangled_ = v;
	}

	void add_row_names(std::vector<W> &&names)
	{
		if (names.size() != base().rows())
			throw std::invalid_argument(
				qs::format("{} Row size mismatch", LOG_HERE));

		row_names_.names_ = std::move(names);
	}

	void add_col_names(std::vector<U> &&names)
	{
		if (names.size() != base().cols())
			throw std::invalid_argument(
				qs::format("{} Col size mismatch", LOG_HERE));

		col_names_.names_ = std::move(names);
	}

	void dbg_print() const
	{
		for (qt::u32 i = 0; i < rows(); i++) {
			for (qt::u32 j = 0; j < cols(); j++) {
				std::cerr << static_cast<int>(get_value(i, j))
					  << " ";
			}
			std::cerr << "\n";
		}
	}
};

template <typename T, typename U, typename W>
struct ov_matrix
    : public matrix_wrapper<ov_matrix<T, U, W>, full_matrix<T>, T, U, W> {
	using Base =
		matrix_wrapper<ov_matrix, meza::view::full::full<T>, T, U, W>;

	ov_matrix(qt::u32 I, qt::u32 J, T *ptr = nullptr)
	    : Base(full_matrix<T>{I, J, meza::shared::layout::DenseRowMajor,
				  ptr})
	{}
};

template <typename T, typename U, typename W>
struct ref_matrix
    : public matrix_wrapper<ref_matrix<T, U, W>, rr_matrix<T>, T, U, W> {
	using Base = matrix_wrapper<ref_matrix, rr_matrix<T>, T, U, W>;

	ref_matrix(qt::u32 I, qt::u32 J, T *ptr = nullptr)
	    : Base(rr_matrix<T>{I, J, meza::shared::layout::RepeatedRow, ptr})
	{}
};

}; // namespace meza::view
#endif // MZ_MATRIX_NAMED
