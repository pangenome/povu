#ifndef MZ_MATRIX_POOL_HAP_COMP_HPP
#define MZ_MATRIX_POOL_HAP_COMP_HPP

#include <algorithm>
#include <cstddef>
#include <optional>
#include <set>
#include <string_view>
#include <vector>

#include <quilt/types.hpp> // for qt::u32, qt::u8, qt::op_t

#include "meza/pool/split.hpp"		  // for matrix_pool
#include "meza/pool/split_pool_types.hpp" // for ov_mat_t

namespace meza::pool::hap_comp
{

// -------
// aliases
// -------

using meza::pool::comparison_op;
using meza::pool::matrix_pool;
using qt::u32;
using qt::u8;

struct haps_comp_set {
	std::set<qt::up_t<qt::u32>> reversals;
	std::set<qt::up_t<qt::u32>> matches;
	std::set<qt::up_t<qt::u32>> mismatches;
};

/**
 * haplotype comparison matrix
 *
 *
 * a square matrix
 *
 * stores from diagonal 1 upwards
 */
template <typename T>
struct hap_comp_matrix {
public:
	// ----------------
	// helpers (public)
	// ----------------

	/**
	 * calculates the number of comparisons in the k-th diagonal of the
	 * haplotype comparison matrix. This is based on the fact that the
	 * matrix is upper triangular and has H_ rows, so the length of the
	 * k-th diagonal is H_ - k.
	 */
	[[nodiscard]] qt::u32 k_len(qt::u32 k) const
	{
		return this->H_ - k;
	}

	/**
	 * calculates the offset in the haplotype comparison matrix for
	 * a given k value. This is based on the number of comparisons
	 * that come before the k-th diagonal in the upper triangle of
	 * the matrix.
	 */
	[[nodiscard]]
	qt::u32 k_offset(qt::u32 k) const
	{
		return (k * k_len(1)) - triangle_number(k - 1) - (k_len(k));
	}

	void clear()
	{
		qt::u32 N = data_size_;
		std::fill(xor_data_, xor_data_ + N, T{});
		std::fill(sum_data_, sum_data_ + N, T{});
		std::fill(xor_ps_.begin(), xor_ps_.end(), 0);
		std::fill(sum_ps_.begin(), sum_ps_.end(), 0);
	}

	void copy_xor_data()
	{
		xor_ps_.assign(xor_data_, xor_data_ + data_size_);
	}

	void copy_sum_data()
	{
		sum_ps_.assign(sum_data_, sum_data_ + data_size_);
	}

	qt::u32 *xor_ps_ptr()
	{
		return xor_ps_.data();
	}

	[[nodiscard]] qt::u32 pool_offset() const
	{
		return pool_offset_;
	}

	/**
	 * number of elements in the upper triangle of the matrix (excluding
	 * diagonal) multiplied by the number of elements per row in the
	 * filter matrix. This is used to determine how much space to
	 * allocate for the xor and sum results in the pool.
	 */
	[[nodiscard]] qt::u32 cols() const
	{
		return filter_->base().cols();
	}

	/**
	 * also the no. of haplotypes (H)
	 */
	[[nodiscard]] qt::u32 rows() const
	{
		return H_;
	}

	[[nodiscard]] const u8 *get_xor_data() const
	{
		return xor_data_;
	}

	[[nodiscard]] u8 *get_xor_data_mut()
	{
		return xor_data_;
	}

	[[nodiscard]] const u8 *get_sum_data() const
	{
		return sum_data_;
	}

	[[nodiscard]] u8 *get_sum_data_mut()
	{
		return sum_data_;
	}

	[[nodiscard]] qt::u32 capacity() const
	{
		return this->capacity_;
	}

	[[nodiscard]] qt::u32 size() const
	{
		return data_size_;
	}

	std::set<qt::up_t<qt::u32>> explore_reversals()
	{
		this->copy_sum_data();

		auto is_reversal = [&](qt::u32 start, qt::u32 end) -> bool
		{
			if (start == 0)
				return sum_ps_[end] == 0;

			// start > 0
			return (sum_ps_[end] - sum_ps_[start - 1]) == 0;
		};

		// compute a special prefix sum
		//
		// only sum if 3
		qt::u32 running{};
		for (qt::u32 i{}; i < data_size_; ++i) {
			if (sum_ps_[i] > 2)
				running += sum_ps_[i];

			sum_ps_[i] = running;
		}

		// explore the pairs and check for reversals based on the
		// special prefix sum
		std::set<qt::up_t<qt::u32>> reversals;
		qt::u32 J = this->cols();
		qt::u32 K = this->rows();
		for (qt::u32 k{1}; k < K; k++) {
			qt::u32 k_len = this->k_len(k);
			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
				auto [h1, h2] = comp_hap_pair(k, k_off);

				const auto &b = filter_->base();
				if (b.is_row_blank(h1) || b.is_row_blank(h2))
					continue;

				qt::u32 start = (k_offset(k) + k_off) * J;
				qt::u32 end = (start + J) - 1;

				if (is_reversal(start, end))
					reversals.emplace(h1, h2);
			}
		}

		return reversals;
	}

	[[nodiscard]] haps_comp_set explore_pairs()
	{
		std::set<qt::up_t<u32>> reversals = this->explore_reversals();
		std::set<qt::up_t<qt::u32>> matches;
		std::set<qt::up_t<qt::u32>> mismatches;

		auto no_inc = [&](qt::u32 start, qt::u32 end) -> bool
		{
			if (start == 0)
				return xor_ps_[end] == 0;

			// start > 0
			return (xor_ps_[end] - xor_ps_[start - 1]) == 0;
		};

		qt::u32 J = this->cols();
		qt::u32 K = this->rows();
		for (qt::u32 k{1}; k < K; k++) {
			qt::u32 k_len = this->k_len(k);
			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
				auto [h1, h2] = comp_hap_pair(k, k_off);

				const auto &b = filter_->base();
				if (b.is_row_blank(h1) || b.is_row_blank(h2))
					continue;

				qt::u32 start = (k_offset(k) + k_off) * J;
				qt::u32 end = (start + J) - 1;

				if (no_inc(start, end))
					matches.emplace(h1, h2);
				else
					mismatches.emplace(h1, h2);
			}
		}

		return {reversals, matches, mismatches};
	}

	void set_filter(const ov_mat_t *f, qt::u32 pool_offset)
	{
		filter_ = f;
		this->H_ = filter_->rows();
		data_size_ = expected_size(*filter_);
		pool_offset_ = pool_offset;

		if (data_size_ > capacity_) {
			std::string err = quilt::shim::format(
				"Data size exceeds capacity of haplotype "
				"comparison matrix. Expected {} Actual {}",
				data_size_, capacity_);

			throw std::runtime_error(err);
		}
	}

	// ------------
	// constructors
	// ------------

	// delete default constructor
	hap_comp_matrix() = delete;
	hap_comp_matrix(const hap_comp_matrix &) = delete;
	hap_comp_matrix &operator=(const hap_comp_matrix &) = delete;

	static hap_comp_matrix<T> create(std::size_t capacity)
	{
		return hap_comp_matrix<T>{capacity};
	}

private:
	/* ============ private data members ======================== */

	// a reference to the filter matrix
	const ov_mat_t *filter_ = nullptr;

	// the offset in the pool where the haplotype comparison matrix
	// starts.
	// This is used to calculate the correct offsets for
	// accessing the data in the pool when performing comparisons.
	qt::u32 pool_offset_;

	// number of haplotypes
	// also no of rows in the filter matrix
	// also the number of rows and cols in the haplotype comparison
	// matrix
	//
	// matrix the number of rows and cols in the hap comp matrix
	qt::u32 H_;

	// the total number of elements in the haplotype comparison
	// matrix diagonal 1 upwards
	qt::u32 data_size_;

	qt::u32 capacity_;

	// pointers to the data in the pool for the xor and sum results
	// these are used to store the results of the comparisons for
	// each pair of haplotypes. The data is stored in a flattened
	// format, where the comparisons for each pair of haplotypes are
	// stored contiguously in memory. The offsets for accessing the
	// correct data for each pair of haplotypes are calculated based
	// on the k value and the number of comparisons that come before
	// it in the upper triangle of the matrix.
	T *xor_data_;
	T *sum_data_;

	std::vector<qt::u32> xor_ps_; // prefix sum
	std::vector<qt::u32> sum_ps_; // prefix sum

	/* ================= private helper functions ================== */

	// ---------------------
	// constructor (private)
	// ---------------------

	hap_comp_matrix(std::size_t capacity)
	    : pool_offset_(0), H_(0), data_size_(0), capacity_(capacity),
	      xor_data_(new T[capacity]), sum_data_(new T[capacity])
	{}

	// -----------------
	// helpers (private)
	// -----------------
	/**
	 * calculates the number of elements in the upper triangle of an
	 * n x n matrix excluding the diagonal
	 */
	[[nodiscard]]
	static qt::u32 triangle_number(qt::u32 n)
	{
		return n * (n + 1) / 2;
	}

	/**
	 * total filter size
	 */
	[[nodiscard]]
	qt::u32 filter_size() const
	{
		return filter_->rows() * filter_->rows();
	}

	/**
	 * Ensure H is set before calling this fn
	 *
	 * calculates the expected size of the haplotype comparison
	 * matrix based on the number of haplotypes (H) and the number
	 * of elements per row in the filter matrix
	 *
	 */
	[[nodiscard]] qt::u32 expected_size(const ov_mat_t &f)
	{
		qt::u32 exp_comparisons =
			(k_len(1) * H_) - triangle_number(H_ - 1);
		qt::u32 elements_per_row = f.cols();
		qt::u32 exp_size = exp_comparisons * elements_per_row;

		return exp_size;
	}

	/**
	 * @brief Given k & k-offset, calculate the row & col pair in
	 * the matrix
	 */
	static qt::up_t<qt::u32> comp_hap_pair(qt::u32 k, qt::u32 k_offset)
	{
		return {k_offset, k + k_offset}; // {h1, h2}
	}
};

} // namespace meza::pool::hap_comp
#endif // MZ_MATRIX_POOL_HAP_COMP_HPP
