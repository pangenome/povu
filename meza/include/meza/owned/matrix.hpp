#ifndef MZ_MATRIX_HPP
#define MZ_MATRIX_HPP

#include <cstddef>
#include <functional> // for std::invoke
#include <iostream>
#include <queue>  // for queue
#include <set>	  // for set
#include <vector> // for vector

#include <log/location.hpp> // for LOG_HERE
#include <quilt/shim.hpp>   // for qs::contains, qs::format
#include <quilt/types.hpp>  // for qt::u32, qt::u8, qt::op_t

#include "meza/shared/shared.hpp" // for layout

namespace meza::matrix
{

using layout = meza::shared::layout;

[[nodiscard]]
std::size_t get_idx(qt::u32 i, qt::u32 j, layout l, const qt::u32 I,
		    const qt::u32 J);

template <typename T>
struct dense_matrix2d {
	static constexpr layout LO_ = layout::DenseRowMajor;
	qt::u32 I_;
	qt::u32 J_;

	// optional: track which rows and cols have data for quick filtering
	std::vector<qt::u8> row_has_data_;
	std::vector<qt::u8> col_has_data_;

	std::vector<T> data_; // row-major order

	dense_matrix2d(qt::u32 I, qt::u32 J, T init = T{})
	    : I_(I), J_(J), row_has_data_(I, 0), col_has_data_(J, 0),
	      data_(static_cast<std::size_t>(I) * J, init)
	{
		if (init != T{}) {
			std::fill(row_has_data_.begin(), row_has_data_.end(),
				  1);
			std::fill(col_has_data_.begin(), col_has_data_.end(),
				  1);
		}
	}

	[[nodiscard]]
	qt::u32 rows() const
	{
		return this->I_;
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return this->J_;
	}

	/**
	 * a notify method to mark a cell as non-blank
	 * useful when updated without using the `set` method
	 */
	void mark_non_blank(qt::u32 i, qt::u32 j)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} ({}, {})", LOG_HERE, i, j));

		row_has_data_[i] = 1;
		col_has_data_[j] = 1;
	}

	[[nodiscard]]
	bool is_row_blank(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(qs::format(
				"{} Invalid row index: {}", LOG_HERE, i));

		return row_has_data_[i] == 0;
	}

	template <class F>
	void update(qt::u32 i, qt::u32 j, F &&fn)
	{
		auto idx = get_idx(i, j, LO_, I_, J_);
		auto &cell = data_[idx];
		std::invoke(std::forward<F>(fn), cell);

		// If you treat "non-blank" as "was ever written", mark
		// unconditionally:
		mark_non_blank(i, j);

		// return out; // lets fn return something useful if you want
	}

	void set(qt::u32 i, qt::u32 j, T value)
	{
		std::size_t idx = get_idx(i, j, LO_, this->I_, this->J_);
		data_[idx] = value;

		// If value is non-zero, mark row/col as non-blank
		if (value != T{})
			this->mark_non_blank(i, j);
	}

	[[nodiscard]]
	const T &at(qt::u32 i, qt::u32 j) const
	{
		return data_[get_idx(i, j, LO_, this->I_, this->J_)];
	}

	[[nodiscard]]
	T &at_mut(qt::u32 i, qt::u32 j)
	{
		return data_[get_idx(i, j, LO_, this->I_, this->J_)];
	}

	// get row method
	// return a slice into the matrix
	[[nodiscard]]
	std::vector<T> copy_row(qt::u32 i) const
	{
		if (i >= this->rows()) {
			// std::string err = pv_cmp::format(
			//	"{} Invalid row index: {}", LOG_HERE, i);
			throw std::out_of_range("");
		}

		std::size_t offset = static_cast<std::size_t>(i) * this->cols();
		return std::vector<T>(data_.begin() + offset,
				      data_.begin() + offset + this->cols());
	}

	[[nodiscard]]
	dense_matrix2d<T> copy_shape() const
	{
		return dense_matrix2d<T>(this->rows(), this->cols());
	}

	[[nodiscard]]
	std::vector<T> xor_rows(qt::u32 i1, qt::u32 i2) const
	{
		std::vector<T> result;
		qt::u32 J = this->cols();

		result.reserve(this->cols());
		auto [r1_it, _] = this->get_row_it(i1);
		auto [r2_it, __] = this->get_row_it(i2);

		for (qt::u32 j{}; j < J; j++, r1_it++, r2_it++)
			result.push_back(*r1_it ^ *r2_it);

		return result;
	}

	void prefix_sum(std::vector<T> &result) const
	{
		for (std::size_t i = 1; i < result.size(); i++)
			result[i] += result[i - 1];
	}

	[[nodiscard]]
	std::vector<T> xor_rows_and_prefix_sum(qt::u32 i1, qt::u32 i2) const
	{
		std::vector<T> result = this->xor_rows(i1, i2);
		this->prefix_sum(result);

		return result;
	}

	[[nodiscard]]
	std::pair<typename std::vector<T>::const_iterator,
		  typename std::vector<T>::const_iterator>
	get_row_it(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(qs::format(
				"{} Invalid row index: {}", LOG_HERE, i));

		std::size_t offset = static_cast<std::size_t>(i) * this->cols();
		return std::make_pair(data_.begin() + offset,
				      data_.begin() + offset + this->cols());
	}

	void dbg_print(std::ostream &os) const
	{
		const qt::u32 I = this->rows();
		const qt::u32 J = this->cols();

		for (qt::u32 i = 0; i < I; i++) {
			for (qt::u32 j{}; j < J; j++) {
				os << this->at(i, j);
				if (J > 0 && j < J - 1)
					os << "\t";
			}
			os << "\n";
		}
	}
};

template <typename T>
struct repeated_row_matrix2d {
	static constexpr layout LO_ = layout::RepeatedRow;
	qt::u32 logical_I_;
	qt::u32 J_;

	// optional: track which rows and cols have data for quick filtering
	qt::u8 row_has_data_ = 0;
	std::vector<qt::u8> col_has_data_;

	std::vector<T> row0; // only store one row

	/**
	 * I is a logical not actual row count
	 */
	repeated_row_matrix2d(qt::u32 I, qt::u32 J, T init = T{})
	    : logical_I_(I), J_(J), row0(J, init), col_has_data_(J, 0)
	{}

	[[nodiscard]]
	qt::u32 rows() const
	{
		return this->logical_I_;
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return this->J_;
	}

	[[nodiscard]]
	bool is_row_blank(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(qs::format(
				"{} Invalid row index: {}", LOG_HERE, i));

		return this->row_has_data_;
	}

	void mark_non_blank(qt::u32 i, qt::u32 j)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} ({}, {})", LOG_HERE, i, j));

		row_has_data_ = 1;
		col_has_data_[j] = 1;
	}

	void set(qt::u32 i, qt::u32 j, T value)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} out of bounds access {}{}",
					   LOG_HERE, i, j));

		row0[j] = value;

		if (value != T{})
			this->mark_non_blank(i, j);
	}

	[[nodiscard]]
	const T &at(qt::u32 i, qt::u32 j) const
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} out of bounds access {}{}",
					   LOG_HERE, i, j));

		return row0[j];
	}

	// get row method
	// return a slice into the matrix
	[[nodiscard]]
	std::vector<T> copy_row(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(qs::format(
				"{} Invalid row index: {}", LOG_HERE, i));

		return row0;
	}

	[[nodiscard]] std::pair<typename std::vector<T>::const_iterator,
				typename std::vector<T>::const_iterator>
	get_row_it(qt::u32 i) const
	{
		if (i >= this->rows())
			throw std::out_of_range(qs::format(
				"{} Invalid row index: {}", LOG_HERE, i));

		return std::make_pair(row0.begin(), row0.end());
	}

	void dbg_print(std::ostream &os) const
	{
		const qt::u32 I = this->rows();
		const qt::u32 J = this->cols();

		for (qt::u32 i = 0; i < I; i++) {
			for (qt::u32 j{}; j < J; j++) {
				os << this->at(i, j);
				if (J > 0 && j < J - 1)
					os << "\t";
			}
			os << "\n";
		}
	}
};

template <typename T>
struct symmetric_square_matrix2d {
	static constexpr layout LO_ = layout::LowerSymmetricSquare;

	qt::u32 I_;
	qt::u32 J_;
	std::size_t N;

	std::vector<qt::u8> row_has_data_;
	std::vector<qt::u8> col_has_data_;

	std::vector<T> data_; // only store lower triangle

	// --------------
	// constructor(s)
	// --------------

	explicit symmetric_square_matrix2d(qt::u32 I, qt::u32 J, T init = T{})
	    : I_{I}, J_{J}, row_has_data_(I, 0), col_has_data_(J, 0)
	{
		if (I_ != J_)
			throw std::invalid_argument(
				qs::format("{} For symmetric matrix, I and "
					   "J must be equal "
					   "(got I={}, J={})",
					   LOG_HERE, I_, J_));

		if (init != T{}) {
			std::fill(row_has_data_.begin(), row_has_data_.end(),
				  1);
			std::fill(col_has_data_.begin(), col_has_data_.end(),
				  1);
		}

		this->N = static_cast<std::size_t>(I) * (I + 1) / 2;
		this->data_.assign(N, init);
	}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	qt::u32 rows() const
	{
		return this->I_;
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return this->J_;
	}

	[[nodiscard]]
	const T &at(qt::u32 i, qt::u32 j) const
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} out of bounds access [{},{}]",
					   LOG_HERE, i, j));

		return data_[get_idx(i, j, LO_, this->I_, this->J_)];
	}

	// ----
	// setters
	// ----
	void set(qt::u32 i, qt::u32 j, T value)
	{
		if (i >= this->rows() || j >= this->cols())
			throw std::out_of_range(
				qs::format("{} out of bounds access {}{}",
					   LOG_HERE, i, j));

		data_[get_idx(i, j, LO_, this->I_, this->J_)] = value;
	}
};

template <typename T>
std::vector<T> prefix_sum(typename std::vector<T>::const_iterator start_it,
			  typename std::vector<T>::const_iterator end_it,
			  qt::u32 N)
{
	std::vector<T> result;
	result.reserve(N);
	auto curr_it = start_it;
	result.push_back(*curr_it);

	curr_it = std::next(curr_it);
	for (; curr_it != end_it; curr_it++)
		result.push_back(result.back() + *curr_it);

	return result;
}

template <typename T>
std::vector<qt::op_t<qt::u32>>
find_context(const repeated_row_matrix2d<T> &ref_matrix,
	     const dense_matrix2d<T> &result_matrix, qt::u32 i)
{
	qt::u32 J = ref_matrix.cols();

	auto [ref_it, _] = ref_matrix.get_row_it(i);
	auto [res_it, res_it_end] = result_matrix.get_row_it(i);

	if (J == 1)
		return {};

	std::vector<T> sum_row = prefix_sum<T>(res_it, res_it_end, J);

	std::vector<qt::op_t<qt::u32>> bounds;
	std::queue<qt::u32> q; // a buffer to store similar cols

	qt::u32 j{}, u{}, v{};

	for (; j < J; j++, ref_it++, res_it++) {
		if (*ref_it != 0 && *res_it == 0)
			q.push(j);

		if (q.size() > 1) {
			u = q.front();
			q.pop();
			v = q.front();

			// if (u,v) are not adjacent & contain a
			// difference
			if ((v - u) > 1 && sum_row[v] - sum_row[u] > 0)
				bounds.emplace_back(u, v);
		}
	}

	return bounds;
}

/**
 * row wise xor
 */
template <typename T>
dense_matrix2d<T> vector_xor(const repeated_row_matrix2d<T> &m1,
			     const dense_matrix2d<T> &m2,
			     const std::set<qt::u32> &skip_rows)
{
	dense_matrix2d<T> res_vec = m2.copy_shape();

	qt::u32 I = m1.rows();
	qt::u32 J = m1.cols();

	for (qt::u32 i{}; i < I; i++) {
		auto [m1_s_it, m1_e_it] = m1.get_row_it(i);
		auto [m2_s_it, m2_e_it] = m2.get_row_it(i);

		if (qs::contains(skip_rows, i))
			continue;

		auto m1_it = m1_s_it;
		auto m2_it = m2_s_it;

		for (qt::u32 j{}; j < J; j++, m1_it++, m2_it++) {
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

	void set_value(qt::u32 i, qt::u32 j, T value)
	{
		base_mut().set(i, j, value);
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
};

template <typename T, typename U, typename W>
struct ov_matrix
    : public matrix_wrapper<ov_matrix<T, U, W>, dense_matrix2d<T>, T, U, W> {
	using Base = matrix_wrapper<ov_matrix, dense_matrix2d<T>, T, U, W>;

	ov_matrix(qt::u32 I, qt::u32 J) : Base(dense_matrix2d<T>{I, J})
	{}
};

template <typename T, typename U, typename W>
struct rep_matrix : public matrix_wrapper<rep_matrix<T, U, W>,
					  repeated_row_matrix2d<T>, T, U, W> {
	using Base =
		matrix_wrapper<rep_matrix, repeated_row_matrix2d<T>, T, U, W>;

	rep_matrix(qt::u32 I, qt::u32 J) : Base(repeated_row_matrix2d<T>{I, J})
	{}
};

using path_matrix =
	ov_matrix<std::vector<std::string>, std::string, std::string>;

using at_matrix = ov_matrix<qt::u32, qt::u32, qt::u32>;

using depth_matrix = ov_matrix<qt::u32, qt::u32, qt::u32>;

using ref_matrix = rep_matrix<qt::u32, qt::u32, qt::u32>;

}; // namespace meza::matrix
#endif // MZ_MATRIX_HPP
