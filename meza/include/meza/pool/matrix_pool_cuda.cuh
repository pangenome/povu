#ifndef MEZA_MATRIX_POOL_CUDA_CUH
#define MEZA_MATRIX_POOL_CUDA_CUH

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <vector_types.h>

#include "meza/pool/matrix_pool.hpp"
#include "meza/pool/split_pool_types.hpp"
#include "meza/shared/shared.hpp" // for layout
#include "meza/view/view.hpp"
#include "quilt/types.hpp"

namespace meza::pool
{

template <typename T>
struct matrix_pool_cuda {
public:
	const matrix_pool<T> &base() const noexcept
	{
		return host_pool_;
	}

	matrix_pool<T> &base_mut() noexcept
	{
		return host_pool_;
	}

	void clear()
	{
		host_pool_.clear();
		// TODO: clear device memory as well
	}

	T *ref_ptr_mut() noexcept
	{
		return d_ref_store_;
	}

	const T *filter_ptr() const noexcept
	{
		return d_filter_store_;
	}

	T *filter_ptr_mut() noexcept
	{
		return d_filter_store_;
	}

	T *xor_ptr_mut() noexcept
	{
		return d_xor_store_;
	}

	/* ================== public methods ========================  */

	void copy_to_device(cudaStream_t stream = 0)
	{
		const auto used = host_pool_.used();
		const auto offsets = host_pool_.offsets();
		T *base = host_pool_.host_data();

		copy_h2d(d_ref_store_, base + offsets.ref_start, used.ref_used,
			 stream);
		copy_h2d(d_filter_store_, base + offsets.filter_start,
			 used.filter_used, stream);
		copy_h2d(d_xor_store_, base + offsets.xor_start, used.xor_used,
			 stream);
	}

	void copy_region_to_host(pool_region region, cudaStream_t stream = 0)
	{
		const auto used = host_pool_.used();
		const auto offsets = host_pool_.offsets();
		T *base = host_pool_.host_data();

		switch (region) {
		case pool_region::Reference:
			copy_d2h(base + offsets.ref_start, d_ref_store_,
				 used.ref_used, stream);
			break;
		case pool_region::Filter:
			copy_d2h(base + offsets.filter_start, d_filter_store_,
				 used.filter_used, stream);
			break;
		case pool_region::Xor:
			copy_d2h(base + offsets.xor_start, d_xor_store_,
				 used.xor_used, stream);
			break;
		}
	}

	void sync(cudaStream_t stream = 0)
	{
		auto err = cudaStreamSynchronize(stream);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}

	/* ================== constructors ========================  */

	explicit matrix_pool_cuda(matrix_pool<T> &hp) : host_pool_(hp)
	{
		auto lens = host_pool_.default_lengths();
		alloc_device(&d_ref_store_, lens.ref_len);
		alloc_device(&d_filter_store_, lens.filter_len);
		alloc_device(&d_xor_store_, lens.xor_len);
	}

	~matrix_pool_cuda()
	{
		free_device(d_ref_store_);
		free_device(d_filter_store_);
		free_device(d_xor_store_);
		free_device(d_haps_xor_store_);
		free_device(d_haps_sum_store_);
	}

	matrix_pool_cuda(const matrix_pool_cuda &) = delete;
	matrix_pool_cuda &operator=(const matrix_pool_cuda &) = delete;
	matrix_pool_cuda(matrix_pool_cuda &&) = delete;
	matrix_pool_cuda &operator=(matrix_pool_cuda &&) = delete;

private:
	/* ============= private data members ======================== */
	matrix_pool<T> &host_pool_;

	T *d_ref_store_ = nullptr;
	T *d_filter_store_ = nullptr;
	T *d_xor_store_ = nullptr;
	T *d_haps_xor_store_ = nullptr;
	T *d_haps_sum_store_ = nullptr;

	/* ============= private methods ======================== */

	static void free_device(T *ptr) noexcept
	{
		if (ptr)
			cudaFree(ptr);
	}

	static void alloc_device(T **ptr, std::size_t count)
	{
		if (count == 0) {
			*ptr = nullptr;
			return;
		}
		auto err = cudaMalloc(reinterpret_cast<void **>(ptr),
				      count * sizeof(T));
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}

	static void copy_h2d(T *dst, const T *src, std::size_t count,
			     cudaStream_t stream)
	{
		if (!count)
			return;
		auto err = cudaMemcpyAsync(dst, src, count * sizeof(T),
					   cudaMemcpyHostToDevice, stream);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}

	static void copy_d2h(T *dst, const T *src, std::size_t count,
			     cudaStream_t stream)
	{
		if (!count)
			return;
		auto err = cudaMemcpyAsync(dst, src, count * sizeof(T),
					   cudaMemcpyDeviceToHost, stream);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}
};

} // namespace meza::pool

#endif // MEZA_MATRIX_POOL_CUDA_CUH
