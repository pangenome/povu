#ifndef MZ_OPS_CUH
#define MZ_OPS_CUH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <quilt/types.hpp>

#include "meza/ops/ops.hpp"

#if MEZA_USE_CUDA
#include <cuda_runtime.h>
#endif

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
void prefix_sum(T *v, std::size_t len)
{
	meza::cpu_ops::prefix_sum_cpu(v, len);
}

} // namespace meza::cuda_ops
#endif // MZ_OPS_CUH
