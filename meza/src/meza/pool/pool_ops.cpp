#include <cstddef>
#include <quilt/types.hpp>

#include "meza/pool/hap_comp.hpp"

#include "meza/pool/split_pool_types.hpp"

#if MEZA_USE_CUDA
#include "meza/ops/ops.cuh"
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/matrix_pool_cuda.cuh"
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
void run_filter(meza::pool::matrix_pool_cuda<qt::u8> &p, qt::u32 N)
{
	const qt::u8 *d_ref = p.ref_ptr_mut();
	const qt::u8 *d_filter = p.filter_ptr_mut();
	qt::u8 *d_xor = p.xor_ptr_mut();

	meza::cuda_ops::cuda_mat_xor(d_ref, d_filter, d_xor, N);
}

void haps_xor(const meza::pool::matrix_pool_cuda<qt::u8> &p,
	      meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
	      qt::u32 len, qt::u32 col_shift, qt::u32 res_shift,
	      cudaStream_t stream = 0)
{
	const qt::u8 *d_f = p.filter_ptr();
	qt::u8 *d_out = cmp_mat.d_haps_xor();
	qt::u32 p_off = cmp_mat.base().pool_offset();

	meza::cuda_ops::cuda_haps_xor(d_f, d_out, p_off, col_shift, res_shift,
				      len, stream);
}

void haps_sum(const meza::pool::matrix_pool_cuda<qt::u8> &p,
	      meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
	      qt::u32 len, qt::u32 col_shift, qt::u32 res_shift,
	      cudaStream_t stream = 0)
{
	const qt::u8 *d_f = p.filter_ptr();
	qt::u8 *d_out = cmp_mat.d_haps_sum();
	qt::u32 p_off = cmp_mat.base().pool_offset();

	meza::cuda_ops::cuda_haps_sum(d_f, d_out, p_off, col_shift, res_shift,
				      len, stream);
}

void run_in_haps(const meza::pool::matrix_pool_cuda<qt::u8> &p,
		 meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
		 meza::pool::comparison_op op)
{
	qt::u32 J = cmp_mat.base().cols();
	qt::u32 K = cmp_mat.base().rows(); // also the no. of haplotypes (H)

	for (qt::u32 k{1}; k < K; k++) {
		qt::u32 col_shift = k * J;
		qt::u32 xor_shift = cmp_mat.base().k_offset(k) * J;

		// number of elements in the k-th diagonal
		qt::u32 len = cmp_mat.base().k_len(k) * J;

		if (op == comparison_op::bitwise_xor)
			haps_xor(p, cmp_mat, len, col_shift, xor_shift);
		else if (op == comparison_op::sum)
			haps_sum(p, cmp_mat, len, col_shift, xor_shift);
	}

	qt::u32 N = cmp_mat.base().size();
	if (op == comparison_op::bitwise_xor)
		cmp_mat.copy_haps_xor_to_host(N);
	else if (op == comparison_op::sum)
		cmp_mat.copy_haps_sum_to_host(N);

	cmp_mat.sync_device(); // TODO: don't do this in prod

	if (op == comparison_op::sum)
		return;

	cmp_mat.base_mut().copy_xor_data();

	meza::cuda_ops::prefix_sum(cmp_mat.base_mut().xor_ps_ptr(), N);
}

haps_comp_set
handle_set(meza::pool::matrix_pool_cuda<qt::u8> &p,
	   // meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8> &cmp_mat,
	   const meza::pool::ov_mat_t &filter_mat, qt::u32 pool_offset)
{
	std::size_t capacity = 512 * 1024 * 1024;
	// hap_comp_matrix<qt::u8> cmp_mat{
	//	filter_mat,  //
	//	pool_offset, //
	//	capacity     //
	// };

	auto cmp_mat =
		meza::pool::hap_comp::hap_comp_matrix<qt::u8>::create(capacity);
	cmp_mat.set_filter(&filter_mat, pool_offset);
	auto cmp_mat_cuda =
		meza::pool::hap_comp::hap_comp_matrix_cuda<qt::u8>{cmp_mat};

	// TODO: parallelise the run_in_haps calls
	run_in_haps(p, cmp_mat_cuda, comparison_op::bitwise_xor);
	run_in_haps(p, cmp_mat_cuda, comparison_op::sum);

	// std::set<qt::up_t<qt::u32>> reversals = comp_mat.find_reversals();
	std::set<qt::up_t<qt::u32>> reversals{};

	auto [matches, mismatches] = cmp_mat.explore_pairs();

	return {reversals, matches, mismatches};
}

} // namespace meza::pool_ops
