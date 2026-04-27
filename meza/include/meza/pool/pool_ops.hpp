#ifndef MEZA_MATRIX_POOL_OPS_HPP
#define MEZA_MATRIX_POOL_OPS_HPP

#include <quilt/types.hpp>

#include "meza/ops/ops.hpp"
#include "meza/pool/hap_comp.hpp"
#include "meza/pool/split.hpp"
#include "meza/pool/split_pool_types.hpp"

#if MEZA_USE_CUDA
#include "meza/ops/ops.cuh"
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/split_cuda.cuh"
#endif

namespace meza::pool_ops
{
using namespace meza::pool::hap_comp; // for haps_comp_set
using qt::u32;
using qt::u8;

template <typename U>
meza::pool::hap_comp::hap_comp_matrix<u8> &base_mat(U &cmp_mat)
{
#if MEZA_USE_CUDA
	return cmp_mat.base_mut();
#else
	return cmp_mat;
#endif
}

void haps_xor_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u8> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift);

void haps_sum_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u8> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift);

#if MEZA_USE_CUDA
void haps_xor_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u8> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift,
		   cudaStream_t stream = 0);

void haps_sum_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u8> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift,
		   cudaStream_t stream = 0);
#endif

/**
 *  T = meza::pool::matrix_pool_cuda<u8>, or
 *  T = meza::pool::matrix_pool<u8>
 *
 *  U = meza::pool::hap_comp::hap_comp_matrix_cuda<u8>, or
 *  U = meza::pool::hap_comp::hap_comp_matrix<u8>
 */
template <typename T, typename U>
void run_in_haps(const T &p, U &f, meza::pool::comparison_op op)
{
	meza::pool::hap_comp::hap_comp_matrix<u8> &cmp_mat = base_mat(f);
	u32 J = cmp_mat.cols();
	u32 K = cmp_mat.rows(); // also the no. of haplotypes (H)

	for (u32 k{1}; k < K; k++) {
		u32 col_shift = k * J;
		u32 xor_shift = cmp_mat.k_offset(k) * J;

		// number of elements in the k-th diagonal
		u32 len = cmp_mat.k_len(k) * J;

#if MEZA_USE_CUDA
		if (op == comparison_op::bitwise_xor)
			haps_xor_cuda(p, f, len, col_shift, xor_shift);
		else if (op == comparison_op::sum)
			haps_sum_cuda(p, f, len, col_shift, xor_shift);
#else
		if (op == comparison_op::bitwise_xor)
			haps_xor_cpu(p, cmp_mat, len, col_shift, xor_shift);
		else if (op == comparison_op::sum)
			haps_sum_cpu(p, cmp_mat, len, col_shift, xor_shift);
#endif
	}

	u32 N = cmp_mat.size();

#if MEZA_USE_CUDA
	if (op == comparison_op::bitwise_xor)
		f.copy_haps_xor_to_host(N);
	else if (op == comparison_op::sum)
		f.copy_haps_sum_to_host(N);

	f.sync_device(); // TODO: don't do this in prod
#endif
	if (op == comparison_op::sum)
		return;

	cmp_mat.copy_xor_data();

#if MEZA_USE_CUDA
	meza::cuda_ops::prefix_sum(cmp_mat.xor_ps_ptr(), N);
#else
	meza::cpu_ops::prefix_sum_cpu(cmp_mat.xor_ps_ptr(), N);
#endif
}

/**
 *  T = meza::pool::matrix_pool_cuda<u8>, or
 *  T = meza::pool::matrix_pool<u8>
 *
 *  U = meza::pool::hap_comp::hap_comp_matrix_cuda<u8>, or
 *  U = meza::pool::hap_comp::hap_comp_matrix<u8>
 *
 * p is the pool containing the reference and filter matrices
 * f is the hap_comp_matrix containing the results of the comparisons
 * (xor and sum)
 *
 * Returns a haps_comp_set containing the results of the comparisons, including
 * the set of reversals and the matches/mismatches for each pair of haplotypes.
 */
template <typename T, typename U>
haps_comp_set handle_set(T &p, U &f)
{
	// TODO: parallelise the run_in_haps calls
	run_in_haps(p, f, comparison_op::bitwise_xor);
	run_in_haps(p, f, comparison_op::sum);

	return base_mat(f).explore_pairs();
}

} // namespace meza::pool_ops

#endif // MEZA_MATRIX_POOL_OPS_HPP
