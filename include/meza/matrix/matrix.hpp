#ifndef MZ_MATRIX_HPP
#define MZ_MATRIX_HPP

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector> // for vector

#include "povu/common/compat.hpp" // for pv_cmp, format
#include "povu/common/core.hpp"

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
template <typename T, typename U>
struct matrix2d {
private:
	pt::u32 I_;
	pt::u32 J_;
	std::vector<T> data_; // row-major order

	std::vector<U> col_names;

	// optional: track which rows and cols have data for quick filtering
	std::vector<pt::u8> row_has_data_;
	std::vector<pt::u8> col_has_data_;

	[[nodiscard]]
	std::size_t get_idx(pt::u32 i, pt::u32 j) const
	{
		if (i >= this->rows() || j >= this->cols()) {
			std::string err = pv_cmp::format(
				"{} out of bounds access {}{}", MODULE, i, j);
			throw std::out_of_range(err);
		}

		return static_cast<std::size_t>(i) * this->cols() + j;
	}

public:
	matrix2d(pt::u32 I, pt::u32 J)
	    : I_{I}, J_{J}, data_(static_cast<std::size_t>(I) * J),
	      row_has_data_(I, 0), col_has_data_(J, 0)
	{}

	matrix2d(pt::u32 I, pt::u32 J, bool with_default_names) : matrix2d(I, J)
	{
		if (!with_default_names)
			return;

		col_names.reserve(J);
		for (pt::u32 j{}; j < J; j++)
			col_names.emplace_back(j);
	}

	[[nodiscard]]
	pt::u32 rows() const
	{
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

	/** it is advisable to use set */
	[[nodiscard]]
	T &at(pt::u32 i, pt::u32 j)
	{
		return data_[get_idx(i, j)];
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

		row_has_data_[i] = 1;
		col_has_data_[j] = 1;
	}

	void set(pt::u32 i, pt::u32 j, T value)
	{
		std::size_t idx = get_idx(i, j);
		data_[idx] = value;

		// If value is non-zero, mark row/col as non-blank
		if (value != T{}) {
			row_has_data_[i] = 1;
			col_has_data_[j] = 1;
		}
	}

	void add_col_names(std::vector<U> &&names)
	{
		if (names.size() != this->cols()) {
			std::string err = pv_cmp::format(
				"{} Column names size {} does "
				"not match column count {}",
				MODULE, names.size(), this->cols());
			throw std::invalid_argument(err);
		}

		col_names = std::move(names);
	}

	[[nodiscard]]
	const std::vector<U> &get_col_names() const
	{
		return col_names;
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

		std::size_t offset = static_cast<std::size_t>(i) * this->cols();
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
		if (j >= this->cols()) {
			throw std::out_of_range("Invalid column index");
		}

		std::vector<T> col(this->rows());
		for (pt::u32 i{}; i < this->rows(); i++)
			col[i] = data_[get_idx(i, j)];

		return col;
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

		std::size_t offset = static_cast<std::size_t>(i) * this->cols();
		return std::make_pair(data_.begin() + offset,
				      data_.begin() + offset + this->cols());
	}

	void dbg_print(std::ostream &os, bool inc_headers = false) const
	{
		if (inc_headers) {
			os << "\t"; // top-left cell is empty
			for (pt::u32 j{}; j < this->cols(); j++) {
				// Safe check: Use the name if it exists,
				// otherwise use the index
				if (j < col_names.size() &&
				    !col_names[j].empty()) {
					os << col_names[j];
				}
				else {
					os << j;
				}
				os << "\t";
			}
			os << "\n";
		}

		for (pt::u32 i = 0; i < this->rows(); i++) {
			if (inc_headers) // print row header
				os << i << "\t";

			for (pt::u32 j = 0; j < this->cols(); j++) {
				os << at(i, j);
				if (this->cols() > 0 && j < this->cols() - 1)
					os << "\t";
			}
			os << "\n";
		}
	}
};

}; // namespace meza::matrix

#endif // MZ_MATRIX_HPP
