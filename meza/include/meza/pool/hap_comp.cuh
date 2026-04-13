#ifndef MZ_MATRIX_POOL_HAP_COMP_CUH
#define MZ_MATRIX_POOL_HAP_COMP_CUH

#include <set>
#include <string_view>
#include <vector>

#include <quilt/types.hpp> // for qt::u32, qt::u8, qt::op_t

// cuda includes
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>

#include "meza/pool/hap_comp.hpp"	  // for hap_comp
#include "meza/pool/split.hpp"		  // for matrix_pool
#include "meza/pool/split_pool_types.hpp" // for ov_mat_t

namespace meza::pool::hap_comp
{

// -------
// aliases
// -------
template <typename T>
struct hap_comp_matrix_cuda {
public:
	void copy_haps_sum_to_host(std::size_t N) noexcept
	{
		T *dst = this->host_hap_comp_.get_sum_data_mut();
		T *src = this->d_haps_sum_store;
		copy_d2h(dst, src, N);
	}

	void copy_haps_xor_to_host(std::size_t N) noexcept
	{
		T *dst = this->host_hap_comp_.get_xor_data_mut();
		T *src = this->d_haps_xor_store;
		copy_d2h(dst, src, N);
	}

	void copy_haps_sum_to_device() noexcept
	{
		T *src = this->host_hap_comp_.get_sum_data_mut();
		T *dst = this->d_haps_sum_store;
		qt::u32 N = this->base().size();
		copy_h2d(dst, src, N);
	}

	void copy_haps_xor_to_device() noexcept
	{
		T *src = this->host_hap_comp_.get_xor_data_mut();
		T *dst = this->d_haps_xor_store;
		qt::u32 N = this->base().size();
		copy_h2d(dst, src, N);
	}

	void sync_device(cudaStream_t stream = 0)
	{
		cudaStreamSynchronize(stream);
	}

	T *d_haps_xor() const noexcept
	{
		return this->d_haps_xor_store;
	}

	T *d_haps_sum() const noexcept
	{
		return this->d_haps_sum_store;
	}

	// -----------
	// base access
	// -----------

	const hap_comp_matrix<T> &base() const noexcept
	{
		return this->host_hap_comp_;
	}

	hap_comp_matrix<T> &base_mut() noexcept
	{
		return this->host_hap_comp_;
	}

	// ------------
	// constructors
	// ------------

	// delete default constructor
	hap_comp_matrix_cuda() = delete;
	hap_comp_matrix_cuda(const hap_comp_matrix_cuda &) = delete;
	hap_comp_matrix_cuda &operator=(const hap_comp_matrix_cuda &) = delete;

	hap_comp_matrix_cuda(hap_comp_matrix<T> &host_hap_comp)
	    : host_hap_comp_(host_hap_comp)
	{
		qt::u32 N = this->base().capacity();
		alloc_device(&d_haps_xor_store, N);
		alloc_device(&d_haps_sum_store, N);
	}

	~hap_comp_matrix_cuda()
	{
		free_device(d_haps_xor_store);
		free_device(d_haps_sum_store);
	}

private:
	hap_comp_matrix<T> &host_hap_comp_;
	T *d_haps_xor_store = nullptr;
	T *d_haps_sum_store = nullptr;

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
		auto err = cudaMalloc(ptr, count * sizeof(T));
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}

	static void copy_h2d(T *dst, const T *src, std::size_t count)
	{
		if (count == 0)
			return;
		auto err = cudaMemcpy(dst, src, count * sizeof(T),
				      cudaMemcpyHostToDevice);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}

	static void copy_d2h(T *dst, const T *src, std::size_t count)
	{
		if (count == 0)
			return;
		auto err = cudaMemcpy(dst, src, count * sizeof(T),
				      cudaMemcpyDeviceToHost);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));
	}
};

} // namespace meza::pool::hap_comp
#endif // MZ_MATRIX_POOL_HAP_COMP_CUH
