#ifndef MZ_MATRIX_HPP
#define MZ_MATRIX_HPP

#include <cstddef>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector> // for vector

#include "povu/common/compat.hpp" // for pv_cmp, format
#include "povu/common/core.hpp"	  // for pt

namespace meza::matrix
{
inline constexpr std::string_view MODULE = "meza::matrix";

/**
 * @brief A 2D matrix with fixed dimensions and high-performance sparsity
 * tracking.
 * * @note Dimensions (I, J) are fixed at construction.
 * @warning The `is_row_blank` and `is_col_blank` tracking logic is optimized
 * for write-heavy/append-only workloads. It assumes data is only added or
 * updated. If elements are "deleted" (reset to T{}), the row/column will still
 * be marked as non-blank unless the internal tracking flags are manually
 * updated.
 */
template <typename T>
struct matrix2d {
private:
	pt::u32 I_;
	pt::u32 J_;
	std::vector<T> data_; // row-major order

	// optional: track which rows and cols have data for quick filtering
	std::vector<pt::u8> row_has_data_;
	std::vector<pt::u8> col_has_data_;

	bool repeated_rows = false;
	pt::u32 false_row_count = 0; // use invalud

	[[nodiscard]]
	std::size_t get_idx(pt::u32 i, pt::u32 j) const
	{
		if (i >= this->rows() || j >= this->cols()) {
			std::string err = pv_cmp::format(
				"{} out of bounds access [{},{}]", MODULE, i,
				j);
			throw std::out_of_range(err);
		}

		if (repeated_rows)
			return j;
		else
			return static_cast<std::size_t>(i) * this->cols() + j;
	}

public:
	matrix2d(pt::u32 I, pt::u32 J)
	    : I_{I}, J_{J}, data_(static_cast<std::size_t>(I) * J),
	      row_has_data_(I, 0), col_has_data_(J, 0)
	{}

	static matrix2d create_repeated_row(pt::u32 false_row_count, pt::u32 J)
	{
		auto m = matrix2d(1, J);
		m.repeated_rows = true;
		m.false_row_count = false_row_count;

		return m;
	}

	matrix2d clone_shape() const
	{
		matrix2d<T> m(this->rows(), this->cols());

		return m;
	}

	[[nodiscard]]
	pt::u32 rows() const
	{
		if (repeated_rows)
			return false_row_count;

		return this->I_;
	}

	[[nodiscard]]
	pt::u32 cols() const
	{
		return this->J_;
	}

	[[nodiscard]]
	bool is_row_blank(pt::u32 i) const
	{
		if (i >= this->rows()) {
			std::string err = pv_cmp::format(
				"{} Invalid row index: {}", MODULE, i);
			throw std::out_of_range(err);
		}

		if (this->repeated_rows)
			return row_has_data_[0] == 0;
		else
			return row_has_data_[i] == 0;
	}

	[[nodiscard]]
	bool is_col_blank(pt::u32 j) const
	{
		if (j >= this->cols()) {
			std::string err = pv_cmp::format(
				"{} Invalid column index: {}", MODULE, j);
			throw std::out_of_range(err);
		}

		return col_has_data_[j] == 0;
	}

	[[nodiscard]]
	const T &at(pt::u32 i, pt::u32 j) const
	{
		return data_[get_idx(i, j)];
	}

	/**
	 * a notify method to mark a cell as non-blank
	 * useful when updated without using the `set` method
	 */
	void mark_non_blank(pt::u32 i, pt::u32 j)
	{
		if (i >= this->rows() || j >= this->cols()) {
			std::string err = pv_cmp::format(
				"{} out of bounds access {}{}", MODULE, i, j);
			throw std::out_of_range(err);
		}

		if (repeated_rows)
			row_has_data_[0] = 1;
		else
			row_has_data_[i] = 1;

		col_has_data_[j] = 1;
	}

	/** it is advisable to use set
	 * if used call mark_non_blank
	 */
	[[nodiscard]]
	T &at(pt::u32 i, pt::u32 j)
	{
		return data_[get_idx(i, j)];
	}

	void set(pt::u32 i, pt::u32 j, T value)
	{
		std::size_t idx = get_idx(i, j);
		data_[idx] = value;

		// If value is non-zero, mark row/col as non-blank
		if (value != T{}) {
			this->mark_non_blank(i, j);
		}
	}

	[[nodiscard]]
	const std::vector<T> &data() const
	{
		return data_;
	}

	// get row method
	// return a slice into the matrix
	[[nodiscard]]
	std::vector<T> get_row(pt::u32 i) const
	{
		if (i >= this->rows()) {
			std::string err = pv_cmp::format(
				"{} Invalid row index: {}", MODULE, i);
			throw std::out_of_range(err);
		}

		pt::u32 row_idx = i;
		if (repeated_rows)
			row_idx = 0;

		std::size_t offset =
			static_cast<std::size_t>(row_idx) * this->cols();
		return std::vector<T>(data_.begin() + offset,
				      data_.begin() + offset + this->cols());
	}

	// get col method
	// return a slice into the matrix
	/**
	 * data stored in row major order
	 * get_col is thus significantly slower than get_row due to cache misses
	 */
	[[nodiscard]]
	std::vector<T> get_col(pt::u32 j) const
	{
		if (j >= this->cols())
			throw std::out_of_range("Invalid column index");

		std::vector<T> col(this->rows());
		for (pt::u32 i{}; i < this->rows(); i++)
			col[i] = data_[get_idx(i, j)];

		return col;
	}

	[[nodiscard]]
	std::vector<T> xor_rows(pt::u32 i1, pt::u32 i2) const
	{
		std::vector<T> result;
		pt::u32 J = this->cols();

		result.reserve(this->cols());
		auto [r1_it, _] = this->get_row_it(i1);
		auto [r2_it, __] = this->get_row_it(i2);

		for (pt::u32 j{}; j < J; j++, r1_it++, r2_it++)
			result.push_back(*r1_it ^ *r2_it);

		return result;
	}

	void prefix_sum(std::vector<T> &result) const
	{
		for (std::size_t i = 1; i < result.size(); i++)
			result[i] += result[i - 1];
	}

	[[nodiscard]]
	std::vector<T> xor_rows_and_prefix_sum(pt::u32 i1, pt::u32 i2) const
	{
		std::vector<T> result = this->xor_rows(i1, i2);
		this->prefix_sum(result);

		return result;
	}

	/**
	 * @brief Retrieves a specific row from the matrix as a range of
	 * iterators.
	 *
	 * iterators are invalidated if the matrix data is modified.
	 *
	 * This method returns a pair of iterators (begin and end) representing
	 * the elements of the specified row in the matrix. It avoids copying
	 * data, providing direct access to the internal storage. The iterators
	 * define a half-open range `[begin, end)`, where `begin` is the start
	 * of the row, and `end` points to one past the last element.
	 *
	 * Note: The `end` iterator is not dereferenceable.
	 *
	 * @tparam T The type of elements in the matrix.
	 * @param i The 0-based index of the row to retrieve.
	 * @return A pair of iterators pointing to the beginning and one past
	 * the end of the row.
	 *
	 * @throws std::out_of_range If the row index `i` is out of bounds.
	 */
	[[nodiscard]]
	std::pair<typename std::vector<T>::const_iterator,
		  typename std::vector<T>::const_iterator>
	get_row_it(pt::u32 i) const
	{
		if (i >= this->rows()) {
			std::string err = pv_cmp::format(
				"{} Invalid row index: {}", MODULE, i);
			throw std::out_of_range(err);
		}

		pt::u32 row_idx = i;
		if (repeated_rows)
			row_idx = 0;

		std::size_t offset =
			static_cast<std::size_t>(row_idx) * this->cols();
		return std::make_pair(data_.begin() + offset,
				      data_.begin() + offset + this->cols());
	}

	[[nodiscard]]
	std::vector<pt::op_t<pt::u32>> find_context(pt::u32 i) const
	{
		pt::u32 J = this->cols();

		auto [start_it, _] = this->get_row_it(i);

		if (J == 1)
			return {};

		std::vector<pt::op_t<pt::u32>> bounds;
		std::queue<pt::u32> q; // a buffer to store similar cols
		auto it{start_it};
		pt::u32 j{}, u{}, v{};

		for (; j < J; j++, it++) {
			if (*it == 0)
				q.push(j);

			if (q.size() > 1) {
				u = q.front();
				q.pop();
				v = q.front();

				// if (u,v) are not adjacent & contain a
				// difference
				if ((v - u) > 1)
					bounds.emplace_back(u, v);
			}
		}

		return bounds;
	}
};

template <typename T>
std::vector<pt::op_t<pt::u32>> find_context(const matrix2d<T> &ref_matrix,
					    const matrix2d<T> &result_matrix,
					    pt::u32 i)
{
	// matrix2d<T, U> res_vec = m2.clone_shape();

	// pt::u32 I = ref_context.rows();
	pt::u32 J = ref_matrix.cols();

	// pt::u32 J = this->cols();

	// auto [start_it, _] = this->get_row_it(i);

	auto [ref_it, _] = ref_matrix.get_row_it(i);
	auto [res_it, __] = result_matrix.get_row_it(i);

	if (J == 1)
		return {};

	std::vector<pt::op_t<pt::u32>> bounds;
	std::queue<pt::u32> q; // a buffer to store similar cols

	// auto m1_it = ref_it;
	// auto m2_it = m2_s_it;

	pt::u32 j{}, u{}, v{};

	for (; j < J; j++, ref_it++, res_it++) {
		if (*ref_it != 0 && *res_it == 0)
			q.push(j);

		if (q.size() > 1) {
			u = q.front();
			q.pop();
			v = q.front();

			// if (u,v) are not adjacent & contain a
			// difference
			if ((v - u) > 1)
				bounds.emplace_back(u, v);
		}
	}

	return bounds;
}

/**
 * row wise xor
 */
template <typename T>
matrix2d<T> vector_xor(const matrix2d<T> &m1, const matrix2d<T> &m2,
		       const std::set<pt::u32> &skip_rows)
{
	matrix2d<T> res_vec = m2.clone_shape();

	pt::u32 I = m1.rows();
	pt::u32 J = m1.cols();

	for (pt::u32 i{}; i < I; i++) {
		auto [m1_s_it, m1_e_it] = m1.get_row_it(i);
		auto [m2_s_it, m2_e_it] = m2.get_row_it(i);

		if (pv_cmp::contains(skip_rows, i))
			continue;

		auto m1_it = m1_s_it;
		auto m2_it = m2_s_it;

		for (pt::u32 j{}; j < J; j++, m1_it++, m2_it++) {
			auto v_res = *m1_it ^ *m2_it;

			res_vec.set(i, j, static_cast<T>(v_res));
		}
	}

	return res_vec;
}

template <typename U>
struct col_names {
	std::vector<U> names_;
};

template <typename W>
struct row_names {
	std::vector<W> names_;
};

template <typename T, typename U, typename W>
struct ov_matrix {
private:
	explicit ov_matrix(matrix2d<T> &&m) : base_matrix_(std::move(m))
	{}

public:
	matrix2d<T> base_matrix_;
	col_names<U> col_names_;
	row_names<W> row_names_;

	// useful for overlays
	bool is_tangled_ = false;
	pt::u32 max_depth_ = 0;

	// --------------
	// constructor(s)
	// --------------

	// TODO: make private and make disp matrix work in a different way
	ov_matrix(pt::u32 I, pt::u32 J) : base_matrix_{I, J}
	{}

	static ov_matrix create_full(pt::u32 I, pt::u32 J)
	{
		return ov_matrix(I, J);
	}

	static ov_matrix create_repeated_row(pt::u32 false_row_count, pt::u32 J)
	{
		return ov_matrix(
			matrix2d<T>::create_repeated_row(false_row_count, J));
	}

	static ov_matrix
	create_from_base_matrix(matrix2d<T> &&m,
				const ov_matrix<T, U, W> &filter_matrix)
	{
		auto ov_m = ov_matrix(std::move(m));
		auto cn = filter_matrix.get_col_names();
		auto rn = filter_matrix.get_row_names();
		ov_m.add_col_names(std::move(cn));
		ov_m.add_row_names(std::move(rn));

		return ov_m;
	}

	// ---------------------------------
	// forwading methods to matrix2dbase
	// ---------------------------------

	const matrix2d<T> &base() const
	{
		return base_matrix_;
	}

	matrix2d<T> &base_mut()
	{
		return base_matrix_;
	}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	const std::vector<W> &get_col_names() const
	{
		return this->col_names_.names_;
	}

	[[nodiscard]]
	const std::vector<W> &get_row_names() const
	{
		return this->row_names_.names_;
	}

	[[nodiscard]]
	pt::u32 get_max_depth() const
	{
		return this->max_depth_;
	}

	[[nodiscard]]
	bool is_tangled() const
	{
		return this->is_tangled_;
	}

	// ---------
	// setter(s)
	// ---------

	void set_max_depth(pt::u32 max_depth)
	{
		this->max_depth_ = max_depth;
	}

	void set_value(pt::u32 i, pt::u32 j, T value)
	{
		this->base_mut().set(i, j, value);
	}

	void set_tangled(bool v)
	{
		this->is_tangled_ = v;
	}

	void add_row_names(std::vector<U> &&names)
	{
		const pt::u32 I = this->base().rows();
		if (names.size() != I) {
			std::string err =
				pv_cmp::format("{} Row names size {} does "
					       "not match row count {}",
					       MODULE, names.size(), I);
			throw std::invalid_argument(err);
		}

		this->row_names_.names_ = std::move(names);
	}

	void add_col_names(std::vector<U> &&names)
	{
		const pt::u32 J = this->base().cols();
		if (names.size() != J) {
			std::string err =
				pv_cmp::format("{} Column names size {} does "
					       "not match column count {}",
					       MODULE, names.size(), J);
			throw std::invalid_argument(err);
		}

		this->col_names_.names_ = std::move(names);
	}

	void dbg_print(std::ostream &os, bool inc_headers = false) const
	{
		// if (inc_headers) {
		//	os << "\t"; // top-left cell is empty
		//	for (pt::u32 j{}; j < this->cols(); j++) {
		//		// Safe check: Use the name if it exists,
		//		// otherwise use the index
		//		if (j < col_names.size() && !col_names.empty())
		//			os << col_names[j];
		//		else
		//			os << j;

		//		os << "\t";
		//	}
		//	os << "\n";
		// }

		// for (pt::u32 i = 0; i < this->rows(); i++) {
		//	if (inc_headers && !this->row_names.empty())
		//		os << this->row_names[i] << "\t";
		//	else if (inc_headers) // row names is empty
		//		os << i << "\t";

		//	for (pt::u32 j{}; j < this->cols(); j++) {
		//		os << at(i, j);
		//		if (this->cols() > 0 && j < this->cols() - 1)
		//			os << "\t";
		//	}
		//	os << "\n";
		// }
	}
};

using path_matrix =
	ov_matrix<std::vector<std::string>, std::string, std::string>;

// a depth matrix can have values of any ...
// template <typename T, typename U, typename W>
using at_matrix = ov_matrix<pt::u32, pt::u32, pt::u32>;

// a depth matrix can have values of any ...
// template <typename T, typename U, typename W>
using depth_matrix = ov_matrix<pt::u32, pt::u32, pt::u32>;

}; // namespace meza::matrix

#endif // MZ_MATRIX_HPP
