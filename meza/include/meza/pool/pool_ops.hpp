#ifndef MEZA_MATRIX_POOL_OPS_HPP
#define MEZA_MATRIX_POOL_OPS_HPP

#include <quilt/types.hpp>

#include "meza/pool/hap_comp.hpp"

#include "meza/pool/split_pool_types.hpp"

#if MEZA_USE_CUDA
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/split_cuda.cuh"
#else
#include "meza/ops/ops.hpp"
#endif

namespace meza::pool_ops
{
using namespace meza::pool::hap_comp; // for haps_comp_set

/**
 * Run the filter operation on the device using the provided matrix pool.
 *
 * @param p the matrix pool containing the reference, filter, and xor regions
 * @param N the number of elements to process in the filter operation
 */
void run_filter(meza::pool::matrix_pool_cuda<qt::u8> &p, qt::u32 N);

void haps_xor(const meza::pool::matrix_pool_cuda<qt::u8> &p,
	      meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
	      qt::u32 len, qt::u32 col_shift, qt::u32 res_shift,
	      cudaStream_t stream = 0);

void haps_sum(const meza::pool::matrix_pool_cuda<qt::u8> &p,
	      meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
	      qt::u32 len, qt::u32 col_shift, qt::u32 res_shift,
	      cudaStream_t stream = 0);

void run_in_haps(const meza::pool::matrix_pool_cuda<qt::u8> &p,
		 meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
		 meza::pool::comparison_op op);

haps_comp_set
handle_set(meza::pool::matrix_pool_cuda<qt::u8> &p,
	   meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat_cuda);

} // namespace meza::pool_ops

#endif // MEZA_MATRIX_POOL_OPS_HPP
