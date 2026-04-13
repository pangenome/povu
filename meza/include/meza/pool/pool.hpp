#ifndef MEZA_MATRIX_POOL_HPP
#define MEZA_MATRIX_POOL_HPP

#include "meza/pool/hap_comp.hpp"
#include "meza/pool/joint.hpp"
#include "meza/pool/pool_ops.hpp"
#include "meza/pool/split.hpp"
#include <cstddef>

#if MEZA_USE_CUDA
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/split_cuda.cuh"
#endif

#include <chrono>
#include <iostream>

namespace meza::pool
{
using namespace meza::pool::hap_comp; // for haps_comp_set
using namespace meza::pool::joint;    // for joint_pool

/**
 * T for the matrix pool (reference, filter, xor),
 * S for the joint pool (depth matrix)
 */
template <typename T, typename S>
struct pool {
public:
	meza::pool::hap_comp::haps_comp_set
	handle_set(const meza::pool::ov_mat_t &filter_mat, qt::u32 pool_offset)
	{
		meza::pool::hap_comp::haps_comp_set cmp_set;
#if MEZA_USE_CUDA
		cmp_mat_cuda.base_mut().set_filter(&filter_mat, pool_offset);
		cmp_set =
			meza::pool_ops::handle_set(mat_pool_cuda, cmp_mat_cuda);
#endif
		return cmp_set;
	}

	void run_convolutions(qt::u32 pool_j_offset)
	{
#if MEZA_USE_CUDA
		mat_pool_cuda.copy_to_device();
		meza::pool_ops::run_filter(this->mat_pool_cuda, pool_j_offset);
		mat_pool_cuda.copy_region_to_host(meza::pool::pool_region::Xor);
#endif
	}

	[[nodiscard]] bool is_full() const
	{
#if MEZA_USE_CUDA
		return mat_pool_cuda.base().is_full();
#else
		return mat_pool_cpu.is_full();
#endif
	}

	void clear_split_pool()
	{
#if MEZA_USE_CUDA
		mat_pool_cuda.clear();
#else
		mat_pool_cpu.clear();
#endif
	}

	void reset_depth_matrix()
	{
		joint_pool_cpu.clear();
	}

	template <typename U, typename W>
	[[nodiscard]] meza::view::ov_matrix<T, U, W>
	alloc_ov_matrix(qt::u32 I, qt::u32 J, pool_region region)
	{
		return mat_pool_cpu.template alloc_ov_matrix<U, W>(I, J,
								   region);
	}

	meza::pool::joint::full_view<S> alloc_depth_matrix(qt::u32 I, qt::u32 J)
	{
		return joint_pool_cpu.alloc_full(I, J);
	}

	/* ================= constructor ============== */

	/*
	 * u8
	 * --------
	 * 1 byte per u8 value
	 * 1024*1024 = 1,048,576 u8 values are
	 *
	 * u32
	 * --------
	 * 4 bytes per u32 value
	 * (1024*1024) / 4 = 262,144
	 * 262144 u32 values are ~1M
	 * 1024 values of u32 are 1M
	 * 2,621,440 u32 values are ~10M
	 */
	pool(std::size_t split_pool_size_mb = 512,		//
	     std::size_t hap_comp_elements = 1024 * 1024,	//
	     std ::size_t joint_pool_size = 10ull * 1024 * 1024 // 10 MiB
	     )
	    : mat_pool_cpu(matrix_pool<T>::create_from_megabytes(
		      split_pool_size_mb)),				  //
	      cmp_mat_cpu(hap_comp_matrix<T>::create(hap_comp_elements)), //
	      joint_pool_cpu(joint_pool<S>::init(joint_pool_size))	  //
#if MEZA_USE_CUDA
	      ,
	      mat_pool_cuda(mat_pool_cpu),
	      cmp_mat_cuda(hap_comp_matrix_cuda<T>{cmp_mat_cpu})
#endif
	{}

private:
	meza::pool::matrix_pool<T> mat_pool_cpu;
	meza::pool::hap_comp::hap_comp_matrix<T> cmp_mat_cpu;
	meza::pool::joint::joint_pool<S> joint_pool_cpu;

#if MEZA_USE_CUDA
	meza::pool::matrix_pool_cuda<T> mat_pool_cuda;
	meza::pool::hap_comp::hap_comp_matrix_cuda<T> cmp_mat_cuda;
#endif
};

}; // namespace meza::pool

#endif // MEZA_MATRIX_POOL_HPP
