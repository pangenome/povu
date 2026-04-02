#ifndef MZ_OPS_HPP
#define MZ_OPS_HPP

#include "quilt/types.hpp"
#include <cassert>
#include <cstdlib>

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;

void cpu_mat_xor(const qt::u8 *a, const qt::u8 *b, qt::u8 *c, qt::u32 N);

/**
 * in-place prefix sum (inclusive scan)
 */
template <typename T>
void cpu_prefix_sum(T *v, size_t len)
{
	for (size_t i = 1; i < len; ++i)
		v[i] = v[i] + v[i - 1];
}

} // namespace meza::cpu_ops
#endif // MZ_OPS_HPP
