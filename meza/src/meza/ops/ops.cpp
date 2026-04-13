#include <cassert>

#include <quilt/types.hpp>

#include "meza/pool/hap_comp.hpp"	  // for haps_comp_set
#include "meza/pool/split.hpp"		  // for matrix_pool
#include "meza/pool/split_pool_types.hpp" // for ov_mat_t

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;

void cpu_mat_xor(const u8 *a, const u8 *b, u8 *c, u32 N)
{
	for (u32 i = 0; i < N; i++)
		c[i] = a[i] ^ b[i];
}

void run_filter(meza::pool::matrix_pool<u8> &ov_pool, u32 N)
{
	const u8 *ref_ptr = ov_pool.ref_start_ptr();
	const u8 *filter_ptr = ov_pool.filter_start_ptr();
	u8 *xor_ptr = ov_pool.xor_start_ptr();
	cpu_mat_xor(ref_ptr, filter_ptr, xor_ptr, N);
}

meza::pool::hap_comp::haps_comp_set
handle_set(meza::pool::matrix_pool<qt::u8> &ov_pool_cuda,
	   const meza::pool::ov_mat_t &filter_mat, qt::u32 pool_offset)
{
	return {};
}

} // namespace meza::cpu_ops
