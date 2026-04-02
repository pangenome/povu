#ifndef MZ_MATRIX_POOL_HAP_COMP_HPP
#define MZ_MATRIX_POOL_HAP_COMP_HPP

// #include <cstddef>
// #include <ostream>
// #include <stdexcept>
// #include "meza/shared/shared.hpp" // for layout

#include <set>
#include <string_view>
#include <vector>

#if MEZA_USE_CUDA
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <vector_types.h>

#include "meza/ops/ops.cuh"
#else
#include "meza/ops/ops.hpp"
#endif

#include "meza/pool/split.hpp" // for matrix_pool, comparison_op, pool_region3
#include "meza/view/view.hpp"
#include "quilt/types.hpp"

namespace meza::pool::hap_comp
{
inline constexpr std::string_view MODULE = "meza::pool::hap_comp";

// -------
// aliases
// -------

using meza::pool::split::comparison_op;
using meza::pool::split::matrix_pool;
using qt::u32;
using qt::u8;

using ov_mat_t = meza::view::ov_matrix<u8, std::string, std::string>;

struct haps_comp_set {
	std::set<qt::up_t<qt::u32>> reversals;
	std::set<qt::up_t<qt::u32>> matches;
	std::set<qt::up_t<qt::u32>> mismatches;
};

haps_comp_set handle_set(matrix_pool<qt::u8> &ov_pool,
			 const ov_mat_t &filter_mat, qt::u32 pool_offset);

/**
 * haplotype comparison matrix
 *
 *
 * a square matrix
 *
 * stores from diagonal 1 upwards
 */
struct hap_comp_matrix {
private:
	// a reference to the filter matrix
	const ov_mat_t &filter_;

	// the offset in the pool where the haplotype comparison matrix
	// starts.
	// This is used to calculate the correct offsets for
	// accessing the data in the pool when performing comparisons.
	qt::u32 pool_offset_;

	// number of haplotypes
	// also no of rows in the filter matrix
	// also the number of rows and cols in the haplotype comparison matrix
	//
	// matrix the number of rows and cols in the hap comp matrix
	qt::u32 H_;

	// the total number of elements in the haplotype comparison matrix
	// diagonal 1 upwards
	qt::u32 data_size_;

	// pointers to the data in the pool for the xor and sum results
	// these are used to store the results of the comparisons for each
	// pair of haplotypes. The data is stored in a flattened format,
	// where the comparisons for each pair of haplotypes are stored
	// contiguously in memory. The offsets for accessing the correct
	// data for each pair of haplotypes are calculated based on the
	// k value and the number of comparisons that come before it in the
	// upper triangle of the matrix.
	qt::u8 *xor_data_;
	qt::u8 *sum_data_;

	std::vector<qt::u32> xor_ps; // prefix sum

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
		return filter_.rows() * filter_.rows();
	}

	[[nodiscard]]
	qt::u32 k_len(qt::u32 k) const
	{
		return this->H_ - k;
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

public:
	// ------------
	// constructors
	// ------------

	// delete default constructor
	hap_comp_matrix() = delete;

	hap_comp_matrix(qt::u8 *pool_xor_ptr, qt::u8 *pool_sum_ptr,
			const ov_mat_t &filter, u32 pool_offset)
	    : filter_(filter), pool_offset_(pool_offset), H_(filter_.rows()),
	      data_size_(expected_size(filter_)), xor_data_(pool_xor_ptr),
	      sum_data_(pool_sum_ptr)
	{}

	// ----------------
	// helpers (public)
	// ----------------

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

	void run_in_haps(matrix_pool<qt::u8> &ov_pool, comparison_op op,
			 cudaStream_t stream = 0)
	{
		qt::u32 J = filter_.cols();
		qt::u32 K = H_;

		// std::cerr << __func__ << " data size: " << data_size_
		// <<
		// "\n";

		for (qt::u32 k{1}; k < K; k++) {
			qt::u32 col_shift = k * J;
			qt::u32 xor_shift = k_offset(k) * J;
			qt::u32 len = k_len(k) * J;

			if (op == comparison_op::bitwise_xor)
				ov_pool.haps_xor(pool_offset_, len, col_shift,
						 xor_shift, stream);
			else if (op == comparison_op::sum)
				ov_pool.hap_sum(pool_offset_, len, col_shift,
						xor_shift, stream);
		}

		if (op == comparison_op::bitwise_xor)
			ov_pool.copy_haps_xor_to_host(data_size_);
		else if (op == comparison_op::sum)
			ov_pool.copy_haps_sum_to_host(data_size_);

		ov_pool.sync_device(stream); // TODO: don't do this in prod

		if (op == comparison_op::sum)
			return;

		xor_ps.assign(xor_data_, xor_data_ + data_size_);

		meza::cuda_ops::prefix_sum(xor_ps.data(), data_size_);
	}

	qt::op_t<std::set<qt::up_t<qt::u32>>> explore_pairs()
	{
		std::set<qt::up_t<qt::u32>> matches;
		std::set<qt::up_t<qt::u32>> mismatches;

		auto no_inc = [&](qt::u32 start, qt::u32 end) -> bool
		{
			if (start == 0)
				return xor_ps[end] == 0;

			// start > 0
			start--;
			return (xor_ps[end] - xor_ps[start]) == 0;
		};

		qt::u32 J = filter_.cols();
		qt::u32 K = H_;
		for (qt::u32 k{1}; k < K; k++) {
			qt::u32 k_len = this->k_len(k);
			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
				auto [h1, h2] = comp_hap_pair(k, k_off);

				if (filter_.base().is_row_blank(h1) ||
				    filter_.base().is_row_blank(h2)) {
					continue;
				}

				qt::u32 start = (k_offset(k) + k_off) * J;
				qt::u32 end = (start + J) - 1;

				if (no_inc(start, end))
					matches.insert({h1, h2});
				else
					mismatches.insert({h1, h2});
			}
		}

		return {matches, mismatches};
	}

	std::set<qt::up_t<qt::u32>> find_reversals()
	{
		std::set<qt::up_t<qt::u32>> reversals;

		// std::cerr << __func__ << " sum data:\n";
		// for (qt::u32 i{}; i < data_size_; i++)
		//	std::cerr << static_cast<u32>(sum_data_[i]) <<
		//", ";
		// std::cerr << "\n";

		std::vector<qt::u32> mask(data_size_, 0);
		for (qt::u32 i{}; i < data_size_; i++)
			if (sum_data_[i] == 1 || sum_data_[i] == 2)
				mask[i] = sum_data_[i];

		meza::cuda_ops::prefix_sum(mask.data(), data_size_);

		// Precompute the `presence` array
		// - tracks whether `1` or `2` exist up to each index
		std::vector<qt::u32> presence(data_size_, 0);
		for (qt::u32 i{}; i < data_size_; i++)
			presence[i] = (i > 0 ? presence[i - 1] : 0) +
				      (sum_data_[i] == 1 || sum_data_[i] == 2);

		// Lambda for the range check
		auto no_inc = [&mask, &presence](qt::u32 start,
						 qt::u32 end) -> bool
		{
			if (start > 0) {
				// Check if prefix sum is constant in
				// range and `1`/`2` exists
				return (mask[end] - mask[start - 1]) == 0 &&
				       (presence[end] - presence[start - 1]) >
					       0;
			}
			else {
				// Special case when start is 0
				return (mask[end] == 0) && presence[end] > 0;
			}
		};

		// std::cerr << __func__ << " sum data summed:\n";
		// for (qt::u32 i{}; i < data_size_; i++)
		//	std::cerr << static_cast<u32>(mask[i]) << ", ";
		// std::cerr << "\n";

		qt::u32 J = filter_.cols();
		qt::u32 K = H_;
		for (qt::u32 k{1}; k < K; k++) {
			qt::u32 k_len = this->k_len(k);
			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
				auto [h1, h2] = comp_hap_pair(k, k_off);

				if (filter_.base().is_row_blank(h1) ||
				    filter_.base().is_row_blank(h2)) {
					continue;
				}

				qt::u32 start = (k_offset(k) + k_off) * J;
				qt::u32 end = (start + J) - 1;

				if (no_inc(start, end))
					reversals.insert({h1, h2});
			}
		}

		return reversals;
	}
};

} // namespace meza::pool::hap_comp
#endif // MZ_MATRIX_POOL_HAP_COMP_HPP
