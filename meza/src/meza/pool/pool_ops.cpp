#include <quilt/types.hpp>

#include "meza/ops/ops.hpp"
#include "meza/pool/hap_comp.hpp"
#include "meza/pool/pool_ops.hpp"
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

/* ====== ==== */
void haps_xor_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u8> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift)
{
	const u8 *d_f = p.filter_start_ptr_const();
	u8 *d_out = cmp_mat.get_xor_data_mut();
	u32 p_off = cmp_mat.pool_offset();

	meza::cpu_ops::cpu_haps_xor(d_f, d_out, p_off, col_shift, res_shift,
				    len);
}

void haps_sum_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u8> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift)
{
	const u8 *d_f = p.filter_start_ptr_const();
	u8 *d_out = cmp_mat.get_sum_data_mut();
	u32 p_off = cmp_mat.pool_offset();

	meza::cpu_ops::cpu_haps_sum(d_f, d_out, p_off, col_shift, res_shift,
				    len);
}

#if MEZA_USE_CUDA
void haps_xor_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u8> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift, cudaStream_t stream)
{
	const u8 *d_f = p.filter_ptr();
	u8 *d_out = cmp_mat.d_haps_xor();
	u32 p_off = cmp_mat.base().pool_offset();

	meza::cuda_ops::cuda_haps_xor(d_f, d_out, p_off, col_shift, res_shift,
				      len, stream);
}

void haps_sum_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u8> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift, cudaStream_t stream)
{
	const u8 *d_f = p.filter_ptr();
	u8 *d_out = cmp_mat.d_haps_sum();
	u32 p_off = cmp_mat.base().pool_offset();

	meza::cuda_ops::cuda_haps_sum(d_f, d_out, p_off, col_shift, res_shift,
				      len, stream);
}
#endif

} // namespace meza::pool_ops
