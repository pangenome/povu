#ifndef MZ_OPS_CUH
#define MZ_OPS_CUH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#if MEZA_USE_CUDA
#include <cuda_runtime.h>
#endif

#include "quilt/types.hpp"

namespace meza::cuda_ops
{
using qt::u32;
using qt::u8;

void cuda_mat_xor(const qt::u8 *a, const qt::u8 *b, qt::u8 *c, qt::u32 N);

void cuda_haps_xor(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		   u32 res_shift, u32 len, cudaStream_t stream);

void cuda_haps_sum(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		   u32 res_shift, u32 len, cudaStream_t stream);

template <typename T>
void cpu_prefix_sum(T *v, size_t len)
{
	for (size_t i = 1; i < len; ++i)
		v[i] = v[i] + v[i - 1];
}

/**
 * prefix sum in-place on the CPU
 *
 * v: input vector of length len, will be modified in-place to contain the
 * prefix sums
 * v[i] = v[0] + v[1] + ... + v[i]
 * len: length of the input vector
 */
template <typename T>
void prefix_sum(T *v, std::size_t len)
{
	cpu_prefix_sum(v, len);
}

} // namespace meza::cuda_ops
#endif // MZ_OPS_CUH
