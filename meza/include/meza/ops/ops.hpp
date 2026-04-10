#ifndef MZ_OPS_HPP
#define MZ_OPS_HPP

#include <cassert>
#include <cstdlib>

#include <quilt/types.hpp>

#include "meza/pool/hap_comp.hpp"	  // for haps_comp_set
#include "meza/pool/matrix_pool.hpp"	  // for matrix_pool
#include "meza/pool/split_pool_types.hpp" // for ov_mat_t

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;

// void run_filter(meza::pool::matrix_pool<u8> &ov_pool, u32 N);

// meza::pool::hap_comp::haps_comp_set
// handle_set(meza::pool::matrix_pool<qt::u8> &ov_pool_cuda,
//	   const meza::pool::ov_mat_t &filter_mat, qt::u32 pool_offset);

// void cpu_mat_xor(const qt::u8 *a, const qt::u8 *b, qt::u8 *c, qt::u32 N);

/**
 * in-place prefix sum (inclusive scan)
 */
template <typename T>
void cpu_prefix_sum(T *v, size_t len)
{
	for (size_t i = 1; i < len; ++i)
		v[i] = v[i] + v[i - 1];
}

std::set<qt::up_t<qt::u32>> find_reversals();

} // namespace meza::cpu_ops
#endif // MZ_OPS_HPP
