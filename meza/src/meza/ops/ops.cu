#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cuda_runtime.h>

#include "quilt/types.hpp"

namespace meza::cuda_ops
{

using qt::u32;
using qt::u8;

// -------
// kernels
// -------
__global__ void matXor(const u8 *a, const u8 *b, u8 *c, u32 N)
{
	u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	// u32 N = I * J;
	if (tid < N)
		c[tid] = a[tid] ^ b[tid];
}

__global__ void hapsSum(const u8 *f, u8 *c, u32 mat_off, u32 col_shift,
			u32 res_shift, u32 len)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < len) {
		u32 f_off = mat_off + tid;
		c[tid + res_shift] = f[f_off] + f[f_off + col_shift];
	}
}

__global__ void hapsXor(const u8 *f, u8 *c, u32 mat_off, u32 col_shift,
			u32 res_shift, u32 len)
{
	uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < len) {
		u32 f_off = mat_off + tid;
		c[tid + res_shift] = f[f_off] ^ f[f_off + col_shift];
	}
}

// --------------
// host functions
// --------------

void cuda_mat_xor(const u8 *d_a, const u8 *d_b, u8 *d_c, u32 N)
{
	u32 blockSize = 256;
	u32 numBlocks = (N + blockSize - 1) / blockSize;
	matXor<<<numBlocks, blockSize>>>(d_a, d_b, d_c, N);

	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}

void cuda_haps_xor(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		   u32 res_shift, u32 len, cudaStream_t stream)
{
	u32 blockSize = 256;
	u32 numBlocks = (len + blockSize - 1) / blockSize;
	hapsXor<<<numBlocks, blockSize, 0, stream>>>(d_f, d_out, mat_off,
						     col_shift, res_shift, len);
	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}

void cuda_haps_sum(const u8 *d_f, u8 *d_out, u32 mat_off, u32 col_shift,
		   u32 res_shift, u32 len, cudaStream_t stream)
{
	u32 blockSize = 256;
	u32 numBlocks = (len + blockSize - 1) / blockSize;
	hapsSum<<<numBlocks, blockSize, 0, stream>>>(d_f, d_out, mat_off,
						     col_shift, res_shift, len);

	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(e));
}
} // namespace meza::cuda_ops
