#ifndef MZ_OPS_HPP
#define MZ_OPS_HPP

#include <cassert>
#include <cstdlib>
#include <set>

#include <quilt/types.hpp>

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;
using qt::up_t;

void cpu_mat_xor(const u8 *a, const u8 *b, u8 *c, u32 N);

void cpu_haps_sum(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		  u32 res_shift, u32 len);

void cpu_haps_xor(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		  u32 res_shift, u32 len);

/**
 * prefix sum in-place on the CPU
 *
 * v: input vector of length len, will be modified in-place to contain the
 * prefix sums
 * v[i] = v[0] + v[1] + ... + v[i]
 * len: length of the input vector
 */
template <typename T>
void prefix_sum_cpu(T *v, std::size_t len)
{
	for (std::size_t i{1}; i < len; ++i)
		v[i] = v[i] + v[i - 1];
}

std::set<up_t<u32>> find_reversals();

} // namespace meza::cpu_ops
#endif // MZ_OPS_HPP
