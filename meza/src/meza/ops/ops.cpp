#include <cassert>

#include <quilt/types.hpp>

namespace meza::cpu_ops
{
using qt::u32;
using qt::u8;

void cpu_mat_xor(const u8 *a, const u8 *b, u8 *c, u32 N)
{
	for (u32 i{}; i < N; i++)
		c[i] = a[i] ^ b[i];
}

void cpu_haps_sum(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		  u32 res_shift, u32 len)
{
	for (u32 i{}; i < len; i++) {
		u32 s{mat_off + i};
		d_out[res_shift + i] = d_f[s] + d_f[s + col_shift];
	}
}

void cpu_haps_xor(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		  u32 res_shift, u32 len)
{
	for (u32 i{}; i < len; i++) {
		u32 s{mat_off + i};
		d_out[res_shift + i] = d_f[s] ^ d_f[s + col_shift];
	}
}

} // namespace meza::cpu_ops
