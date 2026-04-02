#include "quilt/types.hpp"
#include <cassert>

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;

void cpu_mat_xor(const u8 *a, const u8 *b, u8 *c, u32 N)
{
	for (u32 i = 0; i < N; i++)
		c[i] = a[i] ^ b[i];
}

} // namespace meza::cpu_ops
