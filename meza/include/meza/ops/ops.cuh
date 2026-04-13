#ifndef MZ_OPS_CUH
#define MZ_OPS_CUH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <quilt/types.hpp>

#if MEZA_USE_CUDA
#include <cuda_runtime.h>
#endif

#include "meza/pool/hap_comp.hpp"
#include "meza/pool/split_cuda.cuh"

namespace meza::cuda_ops
{
void cuda_mat_xor(const qt::u8 *d_a, const qt::u8 *d_b, qt::u8 *d_c, qt::u32 N);

void cuda_haps_sum(const qt::u8 *d_f, qt::u8 *d_out, qt::u32 mat_off,
		   qt::u32 col_shift, qt::u32 res_shift, qt::u32 len,
		   cudaStream_t stream);

void cuda_haps_xor(const qt::u8 *d_f, qt::u8 *d_out, qt::u32 mat_off,
		   qt::u32 col_shift, qt::u32 res_shift, qt::u32 len,
		   cudaStream_t stream);

template <typename T>
void cpu_prefix_sum(T *v, std::size_t len)
{
	for (std::size_t i = 1; i < len; ++i)
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
