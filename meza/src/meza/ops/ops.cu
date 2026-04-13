#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cuda_runtime.h>

#include "quilt/types.hpp"

#include "meza/pool/hap_comp.hpp"
#include "meza/pool/split_cuda.cuh"
#include "meza/pool/split_pool_types.hpp"

namespace meza::cuda_ops
{

/*
 *
 */

// -------
// kernels
// -------
__global__ void matXor(const qt::u8 *a, const qt::u8 *b, qt::u8 *c, qt::u32 N)
{
	qt::u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	// u32 N = I * J;
	if (tid < N)
		c[tid] = a[tid] ^ b[tid];
}

// --------------
// host functions
// --------------

void cuda_mat_xor(const qt::u8 *d_a, const qt::u8 *d_b, qt::u8 *d_c, qt::u32 N)
{
	qt::u32 blockSize = 256;
	qt::u32 numBlocks = (N + blockSize - 1) / blockSize;
	matXor<<<numBlocks, blockSize>>>(d_a, d_b, d_c, N);

	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}

/*
 *
 */

// -------
// kernels
// -------

__global__ void hapsSum(const qt::u8 *f, qt::u8 *c, qt::u32 mat_off,
			qt::u32 col_shift, qt::u32 res_shift, qt::u32 len)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < len) {
		qt::u32 f_off = mat_off + tid;
		c[tid + res_shift] = f[f_off] + f[f_off + col_shift];
	}
}

__global__ void hapsXor(const qt::u8 *f, qt::u8 *c, qt::u32 mat_off,
			qt::u32 col_shift, qt::u32 res_shift, qt::u32 len)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < len) {
		qt::u32 f_off = mat_off + tid;
		c[tid + res_shift] = f[f_off] ^ f[f_off + col_shift];
	}
}

// --------------
// host functions
// --------------

void cuda_haps_xor(const qt::u8 *d_f, qt::u8 *d_out, qt::u32 mat_off,
		   qt::u32 col_shift, qt::u32 res_shift, qt::u32 len,
		   cudaStream_t stream)
{
	qt::u32 blockSize = 256;
	qt::u32 numBlocks = (len + blockSize - 1) / blockSize;
	hapsXor<<<numBlocks, blockSize, 0, stream>>>(d_f, d_out, mat_off,
						     col_shift, res_shift, len);
	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}

void cuda_haps_sum(const qt::u8 *d_f, qt::u8 *d_out, qt::u32 mat_off,
		   qt::u32 col_shift, qt::u32 res_shift, qt::u32 len,
		   cudaStream_t stream)
{
	qt::u32 blockSize = 256;
	qt::u32 numBlocks = (len + blockSize - 1) / blockSize;
	hapsSum<<<numBlocks, blockSize, 0, stream>>>(d_f, d_out, mat_off,
						     col_shift, res_shift, len);

	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}

} // namespace meza::cuda_ops
